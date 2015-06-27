#!/bin/bash

set -e

size=1024
max=16 # in this script not actually used as max tasks, but just to compute sizes
dt=1
steps=1000
tasks=32

time=(01:30:00 01:45:00 02:00:00 02:15:00 02:30:00 02:45:00 03:00:00 03:15:00 )

function jobfile {
    jobfilename=mpi_lt_$1.job
    echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=$2\n#SBATCH --time=$3\n\necho \"executes on nodes:\"" > $jobfilename
	echo -e 'echo $SLUM_JOB_NODELIST' >> $jobfilename
	echo -e "\necho \"output (if any follows:)\"\n\nif [ ! -z \$1 ]; then\n\tif [ ! -d \"\$1\" ]; then\n\t\tmkdir \"\$1\"\n\tfi\nfi \n\naprun -B /users/stud06/hpcse-adi/mpi/src/main -N $1 --dt $dt --nsteps $steps --localtranspose --benchmark" >> $jobfilename
}

if [ ! -d mpi ]; then
	mkdir mpi
fi
cd mpi


index=1
N=$(echo "scale=4;$size*sqrt($index)"| bc | xargs printf "%1.0f");
jobfile $N $tasks 01:15:00
echo "sbatch mpi_data_$N.job"
sbatch mpi_lt_$N.job

for (( index = 2; index <= $max; index+=2 ))
do
    N=$(echo "scale=4;$size*sqrt($index)"| bc | xargs printf "%1.0f");
    jobfile $N $tasks ${time[$index/2 -1]}
    echo "sbatch mpi_data_$N.job"
    sbatch mpi_lt_$N.job
done


