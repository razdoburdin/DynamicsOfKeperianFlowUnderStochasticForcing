# DynamicsOfKeperianFlowUnderStochasticForcing
Code was used for calculations in the astrophysical investigations. When paper with the results is published, link will be added here.
Tested on Debian 9 with gcc 6.3.0

## Requirements
To compile and run the code you need to have C++11 compiler with [Boost library](https://www.boost.org/) being installed as well as make utility.

## Compilation
To compile project run 
```
make
```
in terminal.

## Run calculations
Parameters of calculations are placed in [configs/params.cfg](configs/params.cfg]:
  + q   -- shear rate
  + R   -- Reynolds number
  + R_b -- Bulk Reynolds number (set R_b=inf to disable bulk viscousity)
  + Ct  -- Courant constant for numerical integration
  + Nt  -- Number of threads
  
To start calculations run one of the following commands in terminal:
  + 
  ```
  make steadyStateTransition
  ```
  to calculate spectrum of perturbations by integration of dynamic equations for the set of SFHs (figures 1 and 2 in the paper).
  + 
  ```
  make spectra[kx]
  ```
  to calculate steady state spectrum of perturbations (figures 3, 4, 5, 6 and 11).
  + 
  ```
  make solutionMap[kx,ky]
  ```
  to calculate dependence of perturbations energy and momentum flux in steady state on forcing Kx and Ky (figures 7 and 8).
  + 
  ```
  make solutionMap[ky,kz]
  ```
  to calculate dependence of perturbations energy and momentum flux in steady state on forcing Ky and Kz (figures 9 and 10).
  + 
  ```
  make optimal[R]
  ```
  to calculate maximal flux and energy (figure 12).
  + 
  ```
  make integrationTest
  ```
  to calculate evolution of single SFH without forcing (figure B1).
    
## Licence
This project is licensed under the GNU General Public License v2.0 - see the [LICENSE](LICENSE) file for details

