#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>

/* Simulation parameters */
#define GRID_W          300
#define GRID_H          300
#define N_AGENTS        20000
#define N_INIT_I        20
#define N_STEPS         730
#define INFECTION_FORCE 0.5

/* Agent states */
#define S 0   
#define E 1  
#define I 2   
#define R 3   

typedef struct {
    int x, y;           
    int state;          
    int time_in_state;  
    int dE, dI, dR;     
    int id;             
} Agent;

/* ─── PRNG helpers ──── */

/**
 * Returns a uniform random double in [0, 1) using a per-process seed.
 * Using rand_r instead of rand() ensures that each MPI process has its
 * own independent random stream.
 */
static double rand_uniform(unsigned int *seed) {
    return rand_r(seed) / (double)RAND_MAX;
}


static double neg_exp(double mean, unsigned int *seed) {
    double u = rand_uniform(seed);
    if (u <= 0.0) u = 1e-15;   /* Avoid log(0) */
    return -mean * log(1.0 - u);
}


int main(int argc, char *argv[]) {

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* CLI arguments: ./infect_central [n_steps] [reproducible 0|1] */
    int n_steps      = (argc > 1) ? atoi(argv[1]) : N_STEPS;
    int reproducible = (argc > 2) ? atoi(argv[2]) : 1;

    Agent *all_agents = NULL;
    if (rank == 0) {
        all_agents = (Agent *)malloc(N_AGENTS * sizeof(Agent));
        if (!all_agents) {
            fprintf(stderr, "malloc failed: all_agents\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        unsigned int seed0 = 42;
        for (int i = 0; i < N_AGENTS; i++) {
            all_agents[i].id            = i;
            all_agents[i].x             = rand_r(&seed0) % GRID_W;
            all_agents[i].y             = rand_r(&seed0) % GRID_H;
            all_agents[i].state         = (i < N_AGENTS - N_INIT_I) ? S : I;
            all_agents[i].time_in_state = 0;
            all_agents[i].dE = (int)ceil(neg_exp(3.0,   &seed0));
            all_agents[i].dI = (int)ceil(neg_exp(7.0,   &seed0));
            all_agents[i].dR = (int)ceil(neg_exp(365.0, &seed0));
        }
    }

    /* ── Domain decomposition ── */
    /*
     * Agents are distributed evenly across processes.
     * The last process absorbs any remainder so that sum(my_n) == N_AGENTS.
     */
    int local_n   = N_AGENTS / size;
    int remainder = N_AGENTS % size;
    int my_n      = (rank == size - 1) ? local_n + remainder : local_n;

    Agent *local_agents = (Agent *)malloc(my_n * sizeof(Agent));
    if (!local_agents) {
        fprintf(stderr, "malloc failed: local_agents (rank %d)\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (rank == 0) {
        for (int r = 1; r < size; r++) {
            int r_n   = (r == size - 1) ? local_n + remainder : local_n;
            int r_off = r * local_n;
            MPI_Send(&all_agents[r_off], r_n * sizeof(Agent), MPI_BYTE,
                     r, 0, MPI_COMM_WORLD);
        }
        memcpy(local_agents, all_agents, my_n * sizeof(Agent));
        free(all_agents);
        all_agents = NULL;
    } else {
        MPI_Recv(local_agents, my_n * sizeof(Agent), MPI_BYTE,
                 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    unsigned int my_seed = reproducible
        ? (unsigned int)(42 + rank * 100003)
        : (unsigned int)(time(NULL) + rank * 99991);

    
    int *local_inf_pos  = (int *)malloc(N_AGENTS * 2 * sizeof(int));
    int *global_inf_pos = (int *)malloc(N_AGENTS * 2 * sizeof(int));
    int *all_counts     = (int *)malloc(size * sizeof(int));

    if (!local_inf_pos || !global_inf_pos || !all_counts) {
        fprintf(stderr, "malloc failed: communication buffers (rank %d)\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int *inf_grid = (int *)calloc(GRID_W * GRID_H, sizeof(int));
    if (!inf_grid) {
        fprintf(stderr, "calloc failed: inf_grid (rank %d)\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    long *seir_log = NULL;
    if (rank == 0) {
        seir_log = (long *)calloc(n_steps * 4, sizeof(long));
        if (!seir_log) {
            fprintf(stderr, "calloc failed: seir_log\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

     double t_start = MPI_Wtime();

    for (int step = 0; step < n_steps; step++) {

        for (int i = 0; i < my_n; i++) {
            local_agents[i].x = rand_r(&my_seed) % GRID_W;
            local_agents[i].y = rand_r(&my_seed) % GRID_H;
        }

        int local_inf_count = 0;
        for (int i = 0; i < my_n; i++) {
            if (local_agents[i].state == I) {
                local_inf_pos[local_inf_count * 2]     = local_agents[i].x;
                local_inf_pos[local_inf_count * 2 + 1] = local_agents[i].y;
                local_inf_count++;
            }
        }

        MPI_Allgather(&local_inf_count, 1, MPI_INT,
                      all_counts, 1, MPI_INT, MPI_COMM_WORLD);

        int *send_counts = (int *)malloc(size * sizeof(int));
        int *displs      = (int *)malloc(size * sizeof(int));
        if (!send_counts || !displs) {
            fprintf(stderr, "malloc failed: Allgatherv buffers (rank %d)\n", rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        int total_inf = 0;
        displs[0] = 0;
        for (int r = 0; r < size; r++) {
            send_counts[r] = all_counts[r] * 2;
            total_inf     += all_counts[r];
        }
        for (int r = 1; r < size; r++)
            displs[r] = displs[r - 1] + send_counts[r - 1];

                        
        MPI_Allgatherv(local_inf_pos,  local_inf_count * 2, MPI_INT,
                       global_inf_pos, send_counts, displs,
                       MPI_INT, MPI_COMM_WORLD);

        free(send_counts);
        free(displs);

        memset(inf_grid, 0, GRID_W * GRID_H * sizeof(int));
        for (int k = 0; k < total_inf; k++) {
            int ix = global_inf_pos[k * 2];
            int iy = global_inf_pos[k * 2 + 1];
            inf_grid[ix * GRID_H + iy]++;
        }

        for (int i = 0; i < my_n; i++) {
            Agent *a = &local_agents[i];

            if (a->state == S) {
                
                int ni = 0;
                for (int dx = -1; dx <= 1; dx++) {
                    for (int dy = -1; dy <= 1; dy++) {
                        int nx = (a->x + dx + GRID_W) % GRID_W;
                        int ny = (a->y + dy + GRID_H) % GRID_H;
                        ni += inf_grid[nx * GRID_H + ny];
                    }
                }
                if (ni > 0) {
                    /* Infection probability: p = 1 - exp(-lambda * Ni) */
                    double p = 1.0 - exp(-INFECTION_FORCE * ni);
                    if (rand_uniform(&my_seed) < p) {
                        a->state         = E;
                        a->time_in_state = 0;
                    }
                }
            }
            else if (a->state == E) {
                a->time_in_state++;
                if (a->time_in_state >= a->dE) {
                    a->state         = I;
                    a->time_in_state = 0;
                }
            }
            else if (a->state == I) {
                a->time_in_state++;
                if (a->time_in_state >= a->dI) {
                    a->state         = R;
                    a->time_in_state = 0;
                }
            }
            else if (a->state == R) {
                a->time_in_state++;
                if (a->time_in_state >= a->dR) {
                    a->state         = S;
                    a->time_in_state = 0;
                }
            }
        }

        long local_seir[4]  = {0, 0, 0, 0};
        long global_seir[4] = {0, 0, 0, 0};
        for (int i = 0; i < my_n; i++)
            local_seir[local_agents[i].state]++;

        MPI_Reduce(local_seir, global_seir, 4, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            seir_log[step * 4 + 0] = global_seir[0];
            seir_log[step * 4 + 1] = global_seir[1];
            seir_log[step * 4 + 2] = global_seir[2];
            seir_log[step * 4 + 3] = global_seir[3];
        }
    }

    double t_end = MPI_Wtime();

    if (rank == 0) {
        printf("Steps=%d | Processes=%d | Reproducible=%d | Time=%.4f s\n",
               n_steps, size, reproducible, t_end - t_start);

        /* Write SEIR time-series */
        char fname[128];
        snprintf(fname, sizeof(fname), "data/seir_data_p%d%s.csv",
                 size, reproducible ? "_repro" : "_nonrepro");
        FILE *fp = fopen(fname, "w");
        if (fp) {
            fprintf(fp, "step,S,E,I,R\n");
            for (int s = 0; s < n_steps; s++)
                fprintf(fp, "%d,%ld,%ld,%ld,%ld\n",
                        s,
                        seir_log[s * 4 + 0],
                        seir_log[s * 4 + 1],
                        seir_log[s * 4 + 2],
                        seir_log[s * 4 + 3]);
            fclose(fp);
            printf("SEIR data written to %s\n", fname);
        } else {
            fprintf(stderr, "Warning: could not open %s for writing\n", fname);
        }

        /* Append timing entry */
        FILE *ft = fopen("data/seir_timing.csv", "a");
        if (ft) {
            fseek(ft, 0, SEEK_END);
            if (ftell(ft) == 0)
                fprintf(ft, "Processes,Reproducible,Time\n");
            fprintf(ft, "%d,%d,%.4f\n", size, reproducible, t_end - t_start);
            fclose(ft);
        }

        free(seir_log);
    }

    /* ── Cleanup ────────────────────────────────────────────────────────── */
    free(local_agents);
    free(local_inf_pos);
    free(global_inf_pos);
    free(all_counts);
    free(inf_grid);

    MPI_Finalize();
    return 0;
}