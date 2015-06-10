#!/bin/bash

set -e

size=1024
max=16
dt=1
steps=10000

time=(01:30:00 01:45:00 02:00:00 02:15:00 02:30:00 02:45:00 03:00:00 03:15:00 )

function params {
	# write N to new file
	echo -e "$1" > params_$1.dat
	# copy all but first line from common params.dat
	tail -n 7 params.dat >> params_$1.dat
}

function jobfile {
	echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=1\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=$2\n#SBATCH --time=$3\n\necho \"executes on nodes:\"" > omp_$1.job
	echo -e 'echo $SLUM_JOB_NODELIST' >> omp_$1.job
	echo -e "\necho \"output (if any follows:)\"\n\nif [ ! -z \$1 ]; then\n\tif [ ! -d \"\$1\" ]; then\n\t\tmkdir \"\$1\"\n\tfi\nfi \n\nexport OMP_NUM_THREADS=$2\naprun -B /users/stud06/hpcse-adi/omp/main -N $1 --dt $dt --nsteps $steps -t $2" >> omp_$1.job
}

if [ ! -d weak_scaling ]; then
	mkdir weak_scaling
fi
cp params.dat weak_scaling
cd weak_scaling

for (( threads = 2; threads <= $max; threads+=2 ))
do
	N=$(echo "scale=4;$size*sqrt($threads)"| bc | xargs printf "%1.0f");
	params $N
	jobfile $N $threads ${time[$threads/2 -1]}
	echo "sbatch omp_$N.job"
	sbatch omp_$N.job
done
