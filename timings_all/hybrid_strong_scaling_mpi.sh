#!/bin/bash

set -e

size=2048
maxtasks=260
threads=8
dt=1
steps=1000

time=(04:00:00 03:45:00 03:30:00 03:15:00 03:00:00 02:45:00 02:30:00 02:15:00 02:00:00 01:45:00 01:30:00)

function jobfile {
    jobfilename=hybrid_data_$2_$3_$1_strong.job
    if (("$3" > "8")); then
    	echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=$2\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=$3\n#SBATCH --ntasks-per-core=2\n#SBATCH --time=$4\n\necho \"executes on nodes:\"" > $jobfilename
    else
    	echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=$2\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=$3\n#SBATCH --time=$4\n\necho \"executes on nodes:\"" > $jobfilename
    fi
	echo -e 'echo $SLUM_JOB_NODELIST' >> $jobfilename
	echo -e "\necho \"output (if any follows:)\"\n\nif [ ! -z \$1 ]; then\n\tif [ ! -d \"\$1\" ]; then\n\t\tmkdir \"\$1\"\n\tfi\nfi \n\naprun -B /users/stud06/hpcse-adi/hybrid/src/main -N $1 --dt $dt --nsteps $steps -t $3 --benchmark" >> $jobfilename
}

if [ ! -d strong_scaling ]; then
	mkdir strong_scaling
fi
cd strong_scaling

tasks=(1 2 4 8 16 32 64 128 256)
#for (( tasks = 1; tasks <= $maxtasks; tasks*=2 ))
for index in {0..8}
do
    jobfile $size ${tasks[$index]} $threads ${time[$index]}
    echo "sbatch hybrid_data_${tasks[$index]}_${threads}_${size}_strong.job"
    sbatch hybrid_data_${tasks[$index]}_${threads}_${size}_strong.job
done
