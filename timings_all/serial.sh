#!/bin/bash

set -e

size=1024
max=16
dt=1
steps=1000

time=(01:30:00 01:45:00 02:00:00 02:15:00 02:30:00 02:45:00 03:00:00 03:15:00 )

function jobfile {
    jobfilename=serial_vec_$1.job
	echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=1\n#SBATCH --time=$3\n\necho \"executes on nodes:\"" > $jobfilename
	echo -e 'echo $SLUM_JOB_NODELIST' >> $jobfilename
	echo -e "\necho \"output (if any follows:)\"\n\nif [ ! -z \$1 ]; then\n\tif [ ! -d \"\$1\" ]; then\n\t\tmkdir \"\$1\"\n\tfi\nfi \n\naprun -B /users/stud06/hpcse-adi/serial/src_daint/main_vectorize -N $1 --dt $dt --nsteps $steps --benchmark" >> $jobfilename
}

if [ ! -d serial ]; then
	mkdir serial
fi
cd serial

#for (( threads = 2; threads <= $max; threads+=2 ))
#do
#	N=$(echo "scale=4;$size*sqrt($threads)"| bc | xargs printf "%1.0f");
#	jobfile $N $threads ${time[$threads/2 -1]}
#    echo "sbatch serial_vec_$N.job"
#    sbatch serial_vec_$N.job
#done

N=512
threads=1
time=01:15:00
jobfile $N $threads $time
echo "sbatch serial_vec_$N.job"
sbatch serial_vec_$N.job
