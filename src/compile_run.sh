g++ ryser.cpp -O0 -o r0
g++ ryser.cpp -O1 -o r1
g++ ryser.cpp -O2 -o r2
g++ ryser.cpp -O3 -o r3
g++ ryser.cpp -Ofast -o rf

g++ ryserP.cpp -O0 -o p0
g++ ryserP.cpp -O1 -o p1
g++ ryserP.cpp -O2 -o p2
g++ ryserP.cpp -O3 -o p3
g++ ryserP.cpp -Ofast -o pf

# Loop through the numbers
# for i in 5 10 15 20 25
# do
#     # Run the r0 binary
#     echo "Running r0 for n=$i"
#     ./r0 $i

#     # Run the r1 binary
#     echo "Running r1 for n=$i"
#     ./r1 $i

#     # Run the r2 binary
#     echo "Running r2 for n=$i"
#     ./r2 $i

#     # Run the r3 binary
#     echo "Running r3 for n=$i"
#     ./r3 $i

#     # Run the rf binary
#     echo "Running rf for n=$i"
#     ./rf $i

#     # Run the p0 binary
#     echo "Running p0 for n=$i"
#     ./p0 $i

#     # Run the p1 binary
#     echo "Running p1 for n=$i"
#     ./p1 $i

#     # Run the p2 binary
#     echo "Running p2 for n=$i"
#     ./p2 $i

#     # Run the p3 binary
#     echo "Running p3 for n=$i"
#     ./p3 $i

#     # Run the pf binary
#     echo "Running pf for n=$i"
#     ./pf $i
# done


for i in 30 35 40
do
     # Run the r2 binary
    echo "Running r2 for n=$i"
    ./r2 $i

    # Run the r3 binary
    echo "Running r3 for n=$i"
    ./r3 $i

    # Run the rf binary
    echo "Running rf for n=$i"
    ./rf $i

    # Run the p2 binary
    echo "Running p2 for n=$i"
    ./p2 $i

    # Run the p3 binary
    echo "Running p3 for n=$i"
    ./p3 $i

    # Run the pf binary
    echo "Running pf for n=$i"
    ./pf $i

done