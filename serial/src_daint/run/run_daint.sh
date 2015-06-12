#!/bin/bash -l
#SBATCH --job-name=adi-mpi-1024-lt
#SBATCH --time=00:30:00
#SBATCH --ntasks=20
#SBATCH --output=out.log.j%j
#SBATCH --error=out.log.j%j


for tasks in {2,4,6,8,10,12,14,16,18,20}; do
    aprun -n $tasks ../main -N 1024 --localtranspose --benchmark
done
