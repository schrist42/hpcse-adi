#!/bin/bash

set -e

size=4096
maxtasks=20
maxthreads=16
dt=1
steps=10000

time=(04:00:00 03:45:00 03:30:00 03:15:00 03:00:00 02:45:00 02:30:00 02:15:00 02:00:00 01:45:00 01:30:00)

function jobfile {
    if (("$3" > "8")); then
    	echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=$2\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=$3\n#SBATCH --ntasks-per-core=2\n#SBATCH --time=$4\n\necho \"executes on nodes:\"" > hybrid_lt_$2_$3_strong.job
    else
    	echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=$2\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=$3\n#SBATCH --time=$4\n\necho \"executes on nodes:\"" > hybrid_lt_$2_$3_strong.job
    fi
	echo -e 'echo $SLUM_JOB_NODELIST' >> hybrid_lt_$2_$3_strong.job
	echo -e "\necho \"output (if any follows:)\"\n\nif [ ! -z \$1 ]; then\n\tif [ ! -d \"\$1\" ]; then\n\t\tmkdir \"\$1\"\n\tfi\nfi \n\naprun -B /users/stud06/hpcse-adi/hybrid/src/main -N $1 --dt $dt --nsteps $steps -t $3 --localtranspose --benchmark" >> hybrid_lt_$2_$3_strong.job
}

if [ ! -d strong_scaling ]; then
	mkdir strong_scaling
fi
cd strong_scaling

for (( tasks = 2; tasks <= $maxtasks; tasks+=2 ))
do
#    for (( threads = 2; threads <=$maxthreads; threads+=2 ))
#    do
        threads=$maxthreads
#    	N=$(echo "scale=4;$size*sqrt($tasks)"| bc | xargs printf "%1.0f");
	    jobfile $size $tasks $threads ${time[$tasks/2 -1]}
	    echo "sbatch hybrid_lt_${tasks}_${threads}_strong.job"
	    sbatch hybrid_lt_${tasks}_${threads}_strong.job
#    done
done
