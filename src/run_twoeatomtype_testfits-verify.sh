#!/bin/bash

g++ -g test_chemfit.cpp -o test_chemfit
g++ -g test_chemfit-energy_evaluator.cpp -o test_chemfit-energy_evaluator

for i in test-He_50F test-Ar_50F # test-Ar+He_50F
do
    for j in 4 5 # 3 6 9 12 -- tok 15 min
    do
        for k in 1 2 3 4
        do
            # Do the fit
        
            #./test_chemfit ${i}.xyzf 50 $k $j > test_chemfit-cpp.log

            # { time python test_chemfit.py ;} > test_chemfit-py.log 2>&1 
        
            #paste b.txt force.txt > compare.txt 
        
            #xmgrace -geometry 1200x1000 compare.txt &

            # Do the verification
        
            # use param file to predict forces
            ./test_chemfit-energy_evaluator ${i}.xyzf 50 $k $j test-Ar+He_50F.${k}.${j}/params.txt > result.dat
        
            awk '!/#/&&/Force/ {print $2}' result.dat > result-force.dat
            awk '!/#/&&/Energy/{print $2}' result.dat > result-energy.dat

            paste ${i}.1.${j}/force.txt result-force.dat > compare-result.txt
            
                   
            xmgrace -geometry 1200x1000 compare-result.txt &
            
            sleep 1
        
            dir=${i}.${k}.${j}-verify
        
            mkdir $dir
        
            mv test_chemfit-cpp.log A.txt b.txt force.txt params.txt test_chemfit-py.log compare.txt result.dat result-force.dat result-energy.dat compare-result.txt $dir    
        
            echo "Finished $dir"
            
        done

    done
done