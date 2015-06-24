#!/bin/bash

set -e

size=1024
max=16
dt=1
steps=1000

time=(01:30:00 01:45:00 02:00:00 02:15:00 02:30:00 02:45:00 03:00:00 03:15:00 )

function jobfile {
    jobfilename=mpi_$1.job
    echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=$2\n#SBATCH --time=$3\n\necho \"executes on nodes:\"" > $jobfilename
	echo -e 'echo $SLUM_JOB_NODELIST' >> $jobfilename
	echo -e "\necho \"output (if any follows:)\"\n\nif [ ! -z \$1 ]; then\n\tif [ ! -d \"\$1\" ]; then\n\t\tmkdir \"\$1\"\n\tfi\nfi \n\naprun -B /users/stud06/hpcse-adi/mpi/src/main -N $1 --dt $dt --nsteps $steps --benchmark" >> $jobfilename
}

if [ ! -d mpi ]; then
	mkdir mpi
fi
cd mpi

for (( tasks = 2; tasks <= $max; tasks+=2 ))
do
    N=$(echo "scale=4;$size*sqrt($tasks)"| bc | xargs printf "%1.0f");
    jobfile $N 20 ${time[$tasks/2 -1]}
    echo "sbatch mpi_$N.job"
    sbatch mpi_$N.job
done

tasks=1
N=$(echo "scale=4;$size*sqrt($tasks)"| bc | xargs printf "%1.0f");
jobfile $N 20 01:15:00
echo "sbatch mpi_$N.job"
sbatch mpi_$N.job


