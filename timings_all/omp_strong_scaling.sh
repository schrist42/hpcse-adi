#!/bin/bash

set -e

size=1024
max=16
dt=1
steps=10000

time=(03:15:00 03:00:00 02:45:00 02:30:00 02:15:00 02:00:00 01:45:00 01:30:00)

function jobfile {
    if (("$2" > "8")); then
	    echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=1\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=$2\n#SBATCH --ntasks-per-core=2\n#SBATCH --time=$3\n\necho \"executes on nodes:\"" > omp_$2_$1_strong.job
    else
        echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=1\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=$2\n#SBATCH --time=$3\n\necho \"executes on nodes:\"" > omp_$2_$1_strong.job
    fi
	echo -e 'echo $SLUM_JOB_NODELIST' >> omp_$2_$1_strong.job
	echo -e "\necho \"output (if any follows:)\"\n\nif [ ! -z \$1 ]; then\n\tif [ ! -d \"\$1\" ]; then\n\t\tmkdir \"\$1\"\n\tfi\nfi \n\nexport OMP_NUM_THREADS=$2\naprun -B /users/stud06/hpcse-adi/omp/src/main -N $1 --dt $dt --nsteps $steps -t $2 --benchmark" >> omp_$2_$1_strong.job
}

if [ ! -d strong_scaling ]; then
	mkdir strong_scaling
fi
cd strong_scaling

for (( threads = 2; threads <= $max; threads+=2 ))
do
	jobfile $size $threads ${time[$threads/2 -1]}
	echo "sbatch omp_${threads}_${size}_strong.job"
	sbatch omp_${threads}_${size}_strong.job
done
