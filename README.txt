Written By: Dorran Howell
11/21/2016
Numerical Modeling in Fortran
HW 7

This package runs a thermal convection model in 2D. In this exercise, one of the goals was to test the numerical stability of the model code for variety of input values. The run_test_cases.py python script executes the compiled Fortran code for a variety of inputs.



Note on compiling:

Please use the 64bit flag when compiling to ensure the program runs properly. To compile this program, add the following packages in this order:

gfortran -fdefault-real-8 poisson_solver_module.f95 convect_module.f95 convect_main.f95
