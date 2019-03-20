# DynamicsOfKeperianFlowUnderStochasticForcing
Computations play an important role in modern astrophysics. Unfortunately many of the research papers are hard to reproduce and check due to computational code was not made public. For that reason I made my code public in hope that it will be useful for someone.

I wrote that code for astrophysical investigations. When paper with the results is published, link will be added here. You are welcome to use and modify code or it's components for you needs.

Program has been tested on
+ 2 cores laptop with Debian 9 and gcc 6.3.0
+ 40 cores CPU server with Ubuntu 16.04 and gcc 5.4.0

## Requirements
To compile and run the code you need to have C++11 compiler with [Boost library](https://www.boost.org/) being installed as well as make utility.

## Compilation
To compile project run 
```
make
```
in terminal.

## Run calculations
Parameters of calculations are placed in [configs/params.cfg](configs/params.cfg):
  + q   -- shear rate.
  + R   -- Reynolds number.
  + R_b -- Bulk Reynolds number (set R_b=inf to disable bulk viscousity).
  + Ct  -- Courant constant for numerical integration.
  + Nt  -- Number of CPU threads in use.
  
To start calculations run one of the following commands in terminal:
  + for calculation of perturbations spectrum by integration of dynamic equations for the set of SFHs (figures 1 and 2 in the paper)
  ```
  make steadyStateTransition
  ```
  + for calculation perturbations steady state spectrum (figures 3, 4, 5, 6 and 11).
  ```
  make spectra[kx]
  ```
  + for calculation of perturbations energy and momentum flux dependence in steady state on forcing Kx and Ky (figures 7 and 8).
  ```
  make solutionMap[kx,ky]
  ```
  + for calculation of perturbations energy and momentum flux dependence in steady state on forcing Ky and Kz (figures 9 and 10).
  ```
  make solutionMap[ky,kz]
  ```
  + for calculation of perturbations maximal flux and energy in steady state (figure 12).
  ```
  make optimal[R]
  ```
  + for calculation of single SFH evolution without forcing (figure B1).
  ```
  make integrationTest
  ```  
## Licence
This project is licensed under the GNU General Public License v2.0 - see the [LICENSE](LICENSE) file for details

