#!/bin/bash

set -e

size=1024
dt=1
steps=100000

time=10:00:00 #(03:15:00 03:00:00 02:45:00 02:30:00 02:15:00 02:00:00 01:45:00 01:30:00)

label=(alpha beta gamma delta epsilon xi eta theta iota kappa lambda mu);
k=(0.046 0.046 0.055 0.055 0.055 0.061 0.063 0.060 0.060 0.063 0.065 0.065);
F=(0.007 0.020 0.024 0.030 0.019 0.024 0.035 0.042 0.050 0.050 0.040 0.050);

function jobfile {
    jobfilename=$3_F$4_k$5.job
	echo -e "#!/bin/bash -l\n\n#SBATCH --ntasks=1\n#SBATCH --time=$2\n\necho \"executes on nodes:\"" > $jobfilename
	echo -e 'echo $SLUM_JOB_NODELIST' >> $jobfilename
	echo -e "\necho \"output (if any follows:)\" \n\naprun -B /users/stud06/hpcse-adi/serial/src_daint/main_vectorize -N $1 --dt $dt --nsteps $steps --pngname $3 -F $4 -k $5" >> $jobfilename
}

for i in {0..11}
do
	jobfile $size $time ${label[$i]} ${F[$i]} ${k[$i]}
	echo "sbatch ${label[$i]}_F${F[$i]}_k${k[$i]}.job"
	sbatch ${label[$i]}_F${F[$i]}_k${k[$i]}.job
done
