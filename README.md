This program computes and diagonalizes the Hückel Hamiltonian of a 1D polyene according to different types of conjugation. 


# INPUT INFORMATION

The input file must be named "input.inp". 

The input has four different sections:

- The first line is the "bond type", which specifies whether there are one or two different values for the bond length. There are two possible values for this keyword:
   - alternated (two different values of bond length) 
   - equal (all bonds are equal) 

- The second line refers to the system being linear or a ring. The keywords are "linear" or "ring". 


- The third line specifies if all the atoms are equal or not:
    - "eq_atoms" refers to the system being composed of only one element
    - "alt_atoms" refers to the system having alternated atoms types (i.e., N and C)

- The last line must be the total number of atoms in the system

## Input example
```
alternated -1.d0 -0.1d0
ring
eq_atoms
10
```


# HOW TO RUN THE CODE

First, to compile the code one must run "make" on the terminal. If the command "make" is not installed, it is also possible to compile it by running "gfortran -o huckel mat_and_diag.f90" (or the fortran compiler of choice). 
Then, the code can be executed just by running `./huckel` in the same directory as the input file. 

# DEPENDENCIES

- Fortran compiler (i.e. gfortran) 
- LAPACK library 


