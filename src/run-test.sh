#!/bin/bash

# Generate the fit

g++ -g test_chemfit.cpp
./a.out
time python test_chemfit.py
paste b.txt force.txt > compare.txt 

xmgrace -geometry 1200x1000 compare.txt


# use param file to predict forces

g++ -g test_chemfit-energy_evaluator.cpp
./a.out 
./a.out  | awk '!/#/{print}'> result.dat
paste b.txt result.dat > compare-result.txt

# Compare fits

xmgrace -geometry 1200x1000 compare.txt compare-result.txt 

# Verify recovery of training forces

paste force.txt result.dat > compare-force2result.txt
xmgrace -geometry 1200x1000 compare-force2result.txt
