#!/bin/bash

set -e

size=4096
tasks=32
threads=8
dt=1
steps=1000

time=(01:30:00 01:45:00 02:00:00 02:15:00 02:30:00 02:45:00 03:00:00 03:15:00 03:30:00)
size=(1024 1449 2048 2508 2896 3238 3547 3831 4096)

function jobfile {
    jobfilename=hybrid_lt_$2_$3_$1.job
    if (("$3" > "8")); then
    	echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=$2\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=$3\n#SBATCH --ntasks-per-core=2\n#SBATCH --time=$4\n\necho \"executes on nodes:\"" > $jobfilename
    else
    	echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=$2\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=$3\n#SBATCH --time=$4\n\necho \"executes on nodes:\"" > $jobfilename
    fi
	echo -e 'echo $SLUM_JOB_NODELIST' >> $jobfilename
    echo -e "\necho \"output (if any follows:)\"\n\nif [ ! -z \$1 ]; then\n\tif [ ! -d \"\$1\" ]; then\n\t\tmkdir \"\$1\"\n\tfi\nfi \n\naprun -B /users/stud06/hpcse-adi/hybrid/src/main -N $1 --dt $dt --nsteps $steps -t $3 --localtranspose --benchmark" >> $jobfilename
}

if [ ! -d hybrid ]; then
	mkdir hybrid
fi
cd hybrid

#for (( tasks = 1; tasks <= $maxtasks; tasks*=2 ))
for index in {0..8}
do
    jobfile ${size[$index]} $tasks $threads ${time[$index]}
    echo "sbatch hybrid_lt_${tasks}_${threads}_${size[$index]}.job"
    sbatch hybrid_lt_${tasks}_${threads}_${size[$index]}.job
done
