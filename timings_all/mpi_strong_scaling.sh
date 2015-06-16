#!/bin/bash

set -e

size=4096
max=20
dt=1
steps=10000

time=(01:30:00 01:45:00 02:00:00 02:15:00 02:30:00 02:45:00 03:00:00 03:15:00 03:30:00 03:45:00 04:00:00)

function jobfile {
	echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=$2\n#SBATCH --time=$3\n\necho \"executes on nodes:\"" > mpi_lt_$2_strong.job
	echo -e 'echo $SLUM_JOB_NODELIST' >> mpi_lt_$2_strong.job
	echo -e "\necho \"output (if any follows:)\"\n\nif [ ! -z \$1 ]; then\n\tif [ ! -d \"\$1\" ]; then\n\t\tmkdir \"\$1\"\n\tfi\nfi \n\naprun -B /users/stud06/hpcse-adi/mpi/src/main -N $1 --dt $dt --nsteps $steps --localtranspose --benchmark" >> mpi_lt_$2_strong.job
}

if [ ! -d strong_scaling ]; then
	mkdir strong_scaling
fi
cd strong_scaling

for (( tasks = 2; tasks <= $max; tasks+=2 ))
do
#	N=$(echo "scale=4;$size*sqrt($tasks)"| bc | xargs printf "%1.0f");
	jobfile $size $tasks ${time[$tasks/2 -1]}
	echo "sbatch mpi_lt_${tasks}_strong.job"
	sbatch mpi_lt_${tasks}_strong.job
done
