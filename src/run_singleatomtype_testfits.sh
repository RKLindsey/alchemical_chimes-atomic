#!/bin/bash

# This script can be used to generate and verify parameters for a single-element type system


g++ -g test_chemfit.cpp -o test_chemfit
g++ -g test_chemfit-energy_evaluator.cpp -o test_chemfit-energy_evaluator

for i in test-He_50F # test-Ar_50F # Took 4.5 min to run on my laptop
do
    for j in 18 # 3 6 9 12
    do
        
        # Do the fit
        
        ./test_chemfit ${i}.xyzf 50 1 $j > test_chemfit-cpp.log

        { time python test_chemfit.py ;} > test_chemfit-py.log 2>&1 
        
        paste b.txt force.txt > compare.txt 
        
        xmgrace -geometry 1200x1000 compare.txt &

        # Do the verification
        
        # use param file to predict forces
        ./test_chemfit-energy_evaluator ${i}.xyzf 50 1 $j params.txt > result.dat
        
        awk '!/#/&&/Force/ {print $2}' result.dat > result-force.dat
        awk '!/#/&&/Energy/{print $2}' result.dat > result-energy.dat

        paste force.txt result-force.dat > compare-result.txt
        
        xmgrace -geometry 1200x1000 compare-result.txt &
         
        sleep 1
         
        dir=${i}.1.${j}
        
        mkdir $dir
        
        mv test_chemfit-cpp.log A.txt b.txt force.txt params.txt test_chemfit-py.log compare.txt result.dat result-force.dat result-energy.dat compare-result.txt $dir    
        
        echo "Finished $dir"

    done
done