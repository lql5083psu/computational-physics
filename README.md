# computational-physics
C++ implementations of classic algorithms in computational physics
Python implementation of result plot
## Simulation of electromagnetic wave propogation in a waveguide using FDTD(Finite-difference time-domain method) approach
1. FDTD is a numetical analysis technique used for modeling computational electrodynamics 
https://en.wikipedia.org/wiki/Finite-difference_time-domain_method

## Schrage random number generation and independence verification
1. Use current system as seed for random number generator
2. Generate a random number between 0 and 1 with Schrage algorithm
3. Odd-th iteration as coordinate x, and even-th iteration as coordinate y. Each (x, y) pair forms a point in 2D plane and is written into a file
4. Calculate correlation using formula: <img src="https://latex.codecogs.com/svg.latex?\Large&space;C(l)=\frac{<x_{n}x_{n+1}>-<x_n>^2}{<x_n^2>-<x_n>^2}" title="C(l)=\frac{<x_{n}x_{n+1}>-<x_n>^2}{<x_n^2>-<x_n>^2}"/>
5. Verify the multi-dimensional independence following formula: <img src="https://latex.codecogs.com/svg.latex?\Large&space;{\chi}^2=\sum_{ij=1}^{K_0}\frac{(n_{ij}-N/k_0^2)^2}{N/k_0^2}" title="{\chi}^2=\sum_{ij=1}^{K_0}\frac{(n_{ij}-N/k_0^2)^2}{N/k_0^2}"/>
