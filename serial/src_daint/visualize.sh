N=256
dt=1

label=(beta gamma iota kappa);
k=(0.046 0.055 0.060 0.063);
F=(0.020 0.024 0.050 0.050);

for index in {0..3}; do
#	echo -e "$N\n$NT\n$Du\n$Dv\n$dt\n${k[index]}\n${F[index]}\n${label[index]}" > "${label[index]}.txt"
	echo "Diffusion ${label[index]} with k=${k[index]} and F=${F[index]} ..."
#	./diffusion < "${label[index]}.txt"
    ./main -N $N --dt $dt -k ${k[index]} -F ${F[index]} --visualize
done
