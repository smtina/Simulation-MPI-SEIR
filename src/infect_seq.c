#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define GRID_W          300
#define GRID_H          300
#define N_AGENTS        20000
#define N_INIT_I        20   
#define N_STEPS         730  
#define INFECTION_FORCE 0.5

#define S 0
#define E 1
#define I 2
#define R 3

typedef struct {
    int x, y;           
    int state;          
    int time_in_state;  
    int dE, dI, dR;    
} Agent;

float negExp(float mean) {
    float u = rand() / (float)RAND_MAX;
    if (u <= 0.0f) u = 1e-15f;
    return -mean * logf(1.0f - u);
}

int main() {
    srand(42); 

    Agent *agents = malloc(N_AGENTS * sizeof(Agent));
    if (!agents) { fprintf(stderr, "Error: agents allocation failed\n"); return 1; }

    int (*grid)[GRID_H] = calloc(GRID_W, sizeof(int[GRID_H]));
    if (!grid) { fprintf(stderr, "Error: grid allocation failed\n"); free(agents); return 1; }

    long stats[N_STEPS][4] = {0};

    for (int i = 0; i < N_AGENTS; i++) {
        agents[i].x = rand() % GRID_W;
        agents[i].y = rand() % GRID_H;
        agents[i].state = (i < N_AGENTS - N_INIT_I) ? S : I;
        agents[i].time_in_state = 0;
        agents[i].dE = (int)ceil(negExp(3.0f));    
        agents[i].dI = (int)ceil(negExp(7.0f));    
        agents[i].dR = (int)ceil(negExp(365.0f));  
    }

    printf("Starting Sequential Simulation (%d steps, %d agents)...\n", N_STEPS, N_AGENTS);
    
    clock_t start = clock();

    
    for (int step = 0; step < N_STEPS; step++) {

        for (int x = 0; x < GRID_W; x++)
            for (int y = 0; y < GRID_H; y++)
                grid[x][y] = 0;

        for (int i = 0; i < N_AGENTS; i++) {
            agents[i].x = rand() % GRID_W;
            agents[i].y = rand() % GRID_H;
            
            if (agents[i].state == I)
                grid[agents[i].x][agents[i].y]++;
        }

        for (int i = 0; i < N_AGENTS; i++) {
            Agent *a = &agents[i];
            stats[step][a->state]++;
          
            if (a->state == S) {
                int ni = 0;
                for (int dx = -1; dx <= 1; dx++) {
                    for (int dy = -1; dy <= 1; dy++) {
                        int nx = (a->x + dx + GRID_W) % GRID_W;
                        int ny = (a->y + dy + GRID_H) % GRID_H;
                        ni += grid[nx][ny];
                    }
                }
              
                if (ni > 0) {
                    double p = 1.0 - exp(-INFECTION_FORCE * ni);
                    if ((rand() / (double)RAND_MAX) < p) {
                        a->state = E;
                        a->time_in_state = 0;
                    }
                }
            }
            else if (a->state == E) {
                a->time_in_state++;
                if (a->time_in_state >= a->dE) {
                    a->state = I;
                    a->time_in_state = 0;
                }
            }
            else if (a->state == I) {
                a->time_in_state++;
                if (a->time_in_state >= a->dI) {
                    a->state = R;
                    a->time_in_state = 0;
                }
            }
            else if (a->state == R) {
                a->time_in_state++;
                if (a->time_in_state >= a->dR) {
                    a->state = S;
                    a->time_in_state = 0;
                }
            }
        }
    }

    clock_t end = clock();
    double cpu_time = (double)(end - start) / CLOCKS_PER_SEC;
    
    FILE *f = fopen("data/seir_sequential.csv", "w");
    if (!f) { 
        fprintf(stderr, "Error: Could not open CSV file. Make sure 'data' folder exists.\n"); 
    } else {
        fprintf(f, "step,S,E,I,R\n");
        for (int s = 0; s < N_STEPS; s++)
            fprintf(f, "%d,%ld,%ld,%ld,%ld\n",
                    s, stats[s][0], stats[s][1], stats[s][2], stats[s][3]);
        fclose(f);
    }

    printf("Simulation completed in %.4f seconds.\n", cpu_time);
    printf("Results exported to: data/seir_sequential.csv\n");

    free(agents);
    free(grid);
    return 0;
}