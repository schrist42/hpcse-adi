#!/bin/bash

set -e

size=512
max=260
dt=1
steps=1000

time=(04:00:00 03:45:00 03:30:00 03:15:00 03:00:00 02:45:00 02:30:00 02:15:00 02:00:00 01:45:00 01:30:00)

function jobfile {
    jobfilename=mpi_data_$1_$2.job
	echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=$2\n#SBATCH --time=$3\n\necho \"executes on nodes:\"" > $jobfilename
    echo -e 'echo $SLUM_JOB_NODELIST' >> $jobfilename
	echo -e "\necho \"output (if any follows:)\"\n\nif [ ! -z \$1 ]; then\n\tif [ ! -d \"\$1\" ]; then\n\t\tmkdir \"\$1\"\n\tfi\nfi \n\naprun -B /users/stud06/hpcse-adi/mpi/src/main -N $1 --dt $dt --nsteps $steps --benchmark" >> $jobfilename
}

if [ ! -d strong_scaling ]; then
	mkdir strong_scaling
fi
cd strong_scaling

tasks=(1 2 4 8 16 32 64 128 256)
#for (( tasks = 1; tasks <= $max; tasks*=2 ))
for i in {0..8}
do
#	N=$(echo "scale=4;$size*sqrt($tasks)"| bc | xargs printf "%1.0f");
#    time_index=$(echo 'l($tasks)/l(2)' | bc -l)
#    echo $time_index
	jobfile $size ${tasks[$i]} ${time[$i]}
	echo "sbatch mpi_data_${size}_${tasks[$i]}.job"
	sbatch mpi_data_${size}_${tasks[$i]}.job
done
