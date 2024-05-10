g++ ryserP2.cpp -O3 -o rp2 -fopenmp

# Loop through the numbers
for i in 5 10 15 20 25
do
    for j in 4 8 16 32 64 128
    do
        echo "n_threads=$j, n=$i"
        ./rp2 $i $j

    done
done