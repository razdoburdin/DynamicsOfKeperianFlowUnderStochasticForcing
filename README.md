# DynamicsOfKeperianFlowUnderStochasticForcing
Code was used for calculations in the astrophysical investigations. When paper with the results is published, link will be added here.
Tested on Debian 9 with gcc 6.3.0

## Requirements
To compile and run the code you need to have C++11 compiler with [Boost library](https://www.boost.org/) being installed as well as make utility.

## Compilation
To compile project run 'make' in terminal.

## Run calculations
Parameters of calculations are placed in [configs/params.cfg](configs/params.cfg]:
  + q   -- shear rate
  + R   -- Reynolds number
  + R_b -- Bulk Reynolds number (set R_b=inf to disable bulk viscousity)
  + Ct  -- Courant constant for numerical integration
  + Nt  -- Number of threads
  
To start calculations:
  + Run 'make steadyStateTransition' in terminal
    to calculate spectrum of perturbations by integration of dynamic equations for the set of SFHs (figures 1 and 2 in the paper).
  + Run 'make spectra[kx]' in terminal
    to calculate steady state spectrum of perturbations (figures 3, 4, 5 and 6).

## Licence
This project is licensed under the GNU General Public License v2.0 - see the [LICENSE](LICENSE) file for details

