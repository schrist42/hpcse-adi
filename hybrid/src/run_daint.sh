#!/bin/bash -l
#SBATCH --job-name=adi-hybrid
#SBATCH --time=00:30:00
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --output=out.log.j%j
#SBATCH --error=out.log.j%j

# hyperthreading: --ntasks-per-core=2

for mpi in {2,4,6,8}; do # max = ntasks
    for omp in {1,2,4,6,8}; do # max = cpus-per-task
        export OMP_NUM_THREADS=$omp
        aprun -n $mpi -N $SLURM_NTASKS_PER_NODE -d $omp main --nthreads $omp -N 512
    done
done
