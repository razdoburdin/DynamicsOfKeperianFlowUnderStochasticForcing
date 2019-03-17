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
Parameters of calculations are placed in [configs/params.cfg](configs/params.cfg):
  + q   -- shear rate
  + R   -- Reynolds number
  + R_b -- Bulk Reynolds number (set R_b=inf to disable bulk viscousity)
  + Ct  -- Courant constant for numerical integration
  + Nt  -- Number of threads
  
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

