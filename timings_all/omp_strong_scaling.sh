#!/bin/bash

set -e

size=512
max=16
dt=1
steps=1000

time=(06:15:00 06:00:00 05:45:00 05:30:00 05:15:00 05:00:00 04:45:00 04:30:00)

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

jobfile $size 1 06:30:00
echo "sbatch omp_1_${size}_strong.job"
sbatch omp_1_${size}_strong.job

for (( threads = 2; threads <= $max; threads+=2 ))
do
	jobfile $size $threads ${time[$threads/2 -1]}
	echo "sbatch omp_${threads}_${size}_strong.job"
	sbatch omp_${threads}_${size}_strong.job
done

