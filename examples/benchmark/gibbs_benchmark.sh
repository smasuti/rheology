#!/usr/bin/env bash

# benchmark example

export PATH=$PATH:../../build

mpirun -n 2 rheology << EOF | tee out.dat 
# number of parameters to optimize
5
# i  init lower upper
1  5.30E+00  2.00E+00  7.00E+00
2  4.30E+00  3.00E+00  5.00E+00
3  3.00E-01 -2.00E+00  1.00E+00
4  5.60E+00  5.00E+00  5.80E+00
5  3.60E+00  3.00E+00  5.30E+00
# Number of samples 
200000
# Conditionals to be evaluated in each dimension 
1000
# location of config file
./config_benchmark.dat
EOF
