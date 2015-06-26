#!/bin/bash

set -e

size=4096
maxtasks=32
max=16
dt=1
steps=1000

time=(06:15:00 06:00:00 05:45:00 05:30:00 05:15:00 05:00:00 04:45:00 04:30:00)

function jobfile {
    jobfilename="hybrid_$2_$3_strong.job"
    if (("$3" > "8")); then
	    echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=$2\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=$3\n#SBATCH --ntasks-per-core=2\n#SBATCH --time=$4\n\necho \"executes on nodes:\"" > $jobfilename
    else
        echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=$2\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=$3\n#SBATCH --time=$4\n\necho \"executes on nodes:\"" > $jobfilename
    fi
	echo -e 'echo $SLUM_JOB_NODELIST' >> $jobfilename
	echo -e "\necho \"output (if any follows:)\"\n\nif [ ! -z \$1 ]; then\n\tif [ ! -d \"\$1\" ]; then\n\t\tmkdir \"\$1\"\n\tfi\nfi \n\nexport OMP_NUM_THREADS=$2\naprun -B /users/stud06/hpcse-adi/hybrid/src/main -N $1 --dt $dt --nsteps $steps -t $3 --benchmark" >> $jobfilename
}

if [ ! -d strong_scaling ]; then
	mkdir strong_scaling
fi
cd strong_scaling

jobfile $size $maxtasks 1 06:30:00
echo "sbatch hybrid_${maxtasks}_1_strong.job"
sbatch hybrid_${maxtasks}_1_strong.job

for (( threads = 2; threads <= $max; threads+=2 ))
do
	jobfile $size $maxtasks $threads ${time[$threads/2 -1]}
	echo "sbatch hybrid_${maxtasks}_${threads}_strong.job"
	sbatch hybrid_${maxtasks}_${threads}_strong.job
done

