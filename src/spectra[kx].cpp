// Copyright 2019 Dmitry N. Razdoburdin.

/* This file is part of DOKFUSF. DOKFUSF is a program that calculats dynamic of Keplerian flow under stochastic forcing.
DOKFUSF is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later version. DOKFUSF is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details. You should have received a copy of the GNU General Public License along with DOKFUSF; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA */

#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

#include <boost/numeric/odeint/integrate/integrate.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/format.hpp>

#include "include/Parameters.h"
#include "include/LyapunovEquations.h"
#include "include/integrator.h"

namespace ode = boost::numeric::odeint;

int main(int ac, char **av) {
    Parameters data(ac, av);
    data.output();

    const double ky    =  1.0;
    const double kz    =  0.0;
    const double dkx   =  0.1;
//    const double kxMax = -pow(data.q / data.invRe * ky, 1.0 / 3.0) + dkx;
//    const double kxMin = -pow(data.q / data.invRe * ky, 1.0 / 3.0);
    const double kxMax =  40;
    const double kxMin = -100;

    std::stringstream SpName;
    SpName << boost::format("E(kx) SoundForcing R=%.0le R_b=%.1le ky=%.2lf kz=%.1lf kxMax=%.1lf") \
    % (1.0 / data.invRe) % (1.0 / data.invRe_b) % (ky) % (kz) % (kxMax);
    std::ofstream fSp;
    fSp.open(SpName.str());

    WaveVector k(data, kxMin, ky, kz);
//    LyapunovEquationWithFlatForcing eqForcing(data, k);
    LyapunovEquationWith2DVorticalWhiteForcing eqForcing(data, k);
//    LyapunovEquationWith2DSoundWhiteForcing eqForcing(data, k);
    Matrix C = ZeroMatrix(4, 4);

    double t = 0;
    fSp << k.x(t) << "\t" << trace(C) << "\t" << get_flux(C) << "\t" << eqForcing.forsingPower(t) << "\n";

    int nx = 0;
    while (k.x(t) < kxMax) {
        eqForcing.make_step_forward(C, t);
        if (k.x(t) - kxMin > (nx + 1) * dkx) {
            ++nx;
            fSp << k.x(t) << "\t" << trace(C) << "\t" << get_flux(C) << "\t" << eqForcing.forsingPower(t) << "\n";
        }
    }

    LyapunovEquationWithoutForcing eq(data, k);
    bool finished = false;
    int  nSteps   = 0;
    const int nStepsMin = 10;
    while (finished == false) {
        eq.make_step_forward(C, t);
        if (k.x(t) - kxMin > (nx + 1) * dkx) {
            ++nx;
            fSp << k.x(t) << "\t" << trace(C) << "\t" << get_flux(C) << "\n";
        }
        finished = (k(t).x() > 0) && \
                   (nSteps > nStepsMin) && \
                   (trace(C) < 1) && \
                   (trace(C) < (kxMax - kxMin) * eqForcing.forsingPower(t));
        ++nSteps;
    }

    IntegratorOut iOut = integrateOverX(data, kxMin, kxMax, ky, kz);
    fprintf(stdout, "kz    = %lf\n", kz);
    fprintf(stdout, "E     = %lf\n", iOut.Ex);
    fprintf(stdout, "F     = %lf\n", iOut.Ix);
    fprintf(stdout, "Ein   = %lf\n", iOut.EInx);
    fprintf(stdout, "E/Ein = %lf\n", iOut.Ex / iOut.EInx);
    fprintf(stdout, "F/Ein = %lf\n", iOut.Ix / iOut.EInx);

    return 0;
}
