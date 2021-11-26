#!/usr/bin/env bash

# This is script for plotting the best fit model. 

export PATH=$PATH:../../build

mkdir -p output
bestfit << EOF | tee out.dat 
# number of parameters to optimize
5
# model parameters
1.0e5
40e3
3.0
400e3
30e3
# location of config file
./config_benchmark.dat
# data directory
./data
# Steady-state parameters (A, E and n)
# A-pre-exponentail factor (log10(A)), E-activation energy and n-stress exponent
4.48 535e3 3.6
EOF
