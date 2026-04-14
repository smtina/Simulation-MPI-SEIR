# SEIR Epidemiological Simulation (MPI / C)

This project implements an Agent-based **SEIR epidemiological model** in C, with:

-  Sequential version
-  Parallel version using MPI (Allgather)
-  Automatic graph generation with Python

---

## Model Description

20 000 agents evolve on a **300 × 300 toroidal grid** (periodic boundaries) over 730 simulated days.  
Each agent carries an epidemiological state drawn from the **SEIR compartmental model**:

| State | Meaning | Transition |
|-------|---------|------------|
| **S** | Susceptible | → E with probability `p = 1 − exp(−0.5 × Nᵢ)` |
| **E** | Exposed | → I after `dE` days ~ Exp(3) |
| **I** | Infectious | → R after `dI` days ~ Exp(7) |
| **R** | Recovered | → S after `dR` days ~ Exp(365) |

At each time step, agents jump to a random cell anywhere on the grid (global mobility).  
Infection is checked in the **Moore neighbourhood** (8 surrounding cells + own cell).

---

## 📂 Project Structure
├── src/
│   ├── infect_seq.c        # Sequential baseline
│   └── infect_central.c    # MPI parallel version 
├── data/                   # CSV outputs (git-ignored)
├── graph/                  # Generated plots 
├── bin/                    # Compiled binaries 
├── plots.py                # Performance & SEIR curve ploting
├── Makefile
└── run_test_infect.sh      # Full benchmark pipeline
## ⚙️ Requirements

- GCC
- MPI (OpenMPI recommended)
- Python 3 with:
  - pandas
  - matplotlib
  - numpy

Install Python dependencies:

```bash
make
```
This creates: -bin/infect_seq , -bin/infect_central

## Run

Run all benchmarks:
```bash 
make test 
```

This will:

->Run sequential simulation
->Run MPI (1, 2, 4, 8 processes)
->Save results in /data
->Generate graphs in /graph
## 📊 Generate Graphs Only
make plot

📈 Outputs
->CSV Files (data/)
->seir_sequential.csv
->seir_data_pX_repro.csv
->speedup_allgather.csv
->Graphs (graph/)
->SEIR evolution curves
->Speedup plots


## Parallelization (MPI)

The MPI version:

Splits agents across processes
Uses MPI_Allgather to share infected positions
Uses MPI_Reduce for global counts
##  Performance

Speedup is computed as:

Speedup = T(1 process) / T(N processes)





