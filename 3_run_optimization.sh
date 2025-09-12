#!/bin/bash

# 사용할 CPU 코어 수
N_CORES=16

echo ">>> Optimizing Initial structure..."
mpirun -np $N_CORES pw.x -in initial_opt.pwi > initial_opt.pwo

echo ">>> Optimizing Final structure..."
mpirun -np $N_CORES pw.x -in final_opt.pwi > final_opt.pwo

echo ">>> Optimization complete!"
