#!/bin/bash

set -e

size=1024
max=16
dt=1
steps=1000
threads=8

time=(01:30:00 01:45:00 02:00:00 02:15:00 02:30:00 02:45:00 03:00:00 03:15:00 )

function jobfile {
    if (("$2" > "8")); then
	    echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=1\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=$2\n#SBATCH --ntasks-per-core=2\n#SBATCH --time=$3\n\necho \"executes on nodes:\"" > omp_$1.job
    else
        echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=1\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=$2\n#SBATCH --time=$3\n\necho \"executes on nodes:\"" > omp_$1.job
    fi
	echo -e 'echo $SLUM_JOB_NODELIST' >> omp_$1.job
	echo -e "\necho \"output (if any follows:)\"\n\nif [ ! -z \$1 ]; then\n\tif [ ! -d \"\$1\" ]; then\n\t\tmkdir \"\$1\"\n\tfi\nfi \n\nexport OMP_NUM_THREADS=$2\naprun -B /users/stud06/hpcse-adi/omp/src/main -N $1 --dt $dt --nsteps $steps -t $2 --benchmark" >> omp_$1.job
}

if [ ! -d omp ]; then
	mkdir omp
fi
cd omp

index=1
N=$(echo "scale=4;$size*sqrt($index)"| bc | xargs printf "%1.0f");
jobfile $N $threads 01:15:00
echo "sbatch omp_$N.job"
sbatch omp_$N.job

for (( index = 2; index <= $max; index+=2 ))
do
	N=$(echo "scale=4;$size*sqrt($index)"| bc | xargs printf "%1.0f");
	jobfile $N $threads ${time[$index/2 -1]}
	echo "sbatch omp_$N.job"
	sbatch omp_$N.job
done
