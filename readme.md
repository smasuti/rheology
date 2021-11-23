This is repository for fitting the stress vs strain curves using MCMC method. 
Current implementation includes Gibbs sampling. Gibbs sampling works with and 
without MPI parallelism.

There are two flags to control the usage of MPI (USE_MPI) and the method of numerical 
integration (PEC), which can be enabled/disabled in src/include.f90. The code is 
mainly tested with PEC (predictor-corrector method). If you use Runge-Kutta method,
you may come across the error of step size being too small.Currently only nonlinear 
Burgers model is implemeted. The user is free to implement their own model.

In order to compile the code you will need gfortran and mpif90. If you do not 
want to use the MPI then consider replacing the mpif90 with gfortran or suitable
compiler. The code is mainly tested on MacOS, but can be easily compiled on 
Windows if the suitable libraries and compilers are installed. Once successfully
compiler, you can run the following command: 

Go to the examples directory and execute the bash script to get started.

Make sure to write the config file correctly and keep all the data in the same 
directory. Best place to start would be to familiarise with benchmark example.
There is gibbs_benchmark.sh to run the MCMC inversion and also bestfit.sh to generate
the forward best fit for the current model.

If you use the code or method please cite the following paper: 
Masuti, S., & Barbot, S. MCMC inversion of the transient and steady-state creep 
flow law parameters of dunite under dry and wet conditions, EPS 73, 208 (2021). 
https://doi.org/10.1186/s40623-021-01543-9
