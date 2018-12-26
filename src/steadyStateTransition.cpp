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
#include <vector>

#include <boost/numeric/odeint/integrate/integrate.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/format.hpp>

#include "include/Parameters.h"
#include "include/LyapunovEquations.h"

namespace ode = boost::numeric::odeint;

void spectraOut(const std::vector <WaveVector>& k, const std::vector <Matrix>& C, const Parameters& data, double t) {
    std::stringstream SpName;
    SpName << boost::format("NonSteadySpectra R=%.0le R_b=%.1le ky=%.2lf kz=%.1lf t=%.3lf") \
    % (1.0 / data.invRe) % (1.0 / data.invRe_b) % (k[0].y()) % (k[0].z()) % (t);
    std::ofstream fSp;
    fSp.open(SpName.str());

    for (uint i = 0; i < k.size(); ++i) {
        fSp << k[i].x(t) << "\t" << trace(C[i]) << "\n";
    }

    fSp.close();
}

void addPointOfSingleSFH(const WaveVector& k, const Matrix C, const Parameters& data, double t) {
    std::stringstream SpName;
    SpName << boost::format("SteadySpectra R=%.0le R_b=%.1le ky=%.2lf kz=%.1lf") \
    % (1.0 / data.invRe) % (1.0 / data.invRe_b) % (k.y()) % (k.z());
    std::ofstream fSp;
    fSp.open(SpName.str(), std::ios_base::app);

    fSp << k.x(t) << "\t" << trace(C) << "\n";
    fSp.close();
}

int main(int ac, char **av) {
    Parameters data(ac, av);
    data.output();

    const double ky     =  1.0;
    const double kz     =  0.0;
    const double kxFMax =  3.5;
    const double kxFMin = -3.5;

    const double kxMax  =   20.0;
    const double kxMin  =  -50.0;
    const double dkx    =   0.1;
    const double dt     =   1e-3 / (data.q * ky);

    const int nKx       = static_cast <int> ((kxMax - kxMin) / dkx) + 1;
    std::vector <double> kx;
    for (int iKx = 0; iKx < nKx; ++iKx) {
        kx.push_back(kxMax - dkx * iKx);
    }

    int iKFirst = nKx - static_cast <int> ((kxFMin - kxMin) / dkx) - 1;

    std::vector <WaveVector> k;
    std::vector <LyapunovEquationWithFlatForcing> eqForcing;
    std::vector <LyapunovEquationWithoutForcing> eqFree;
    std::vector <Matrix> C;
    std::vector <double> tMax;
    for (uint i = 0; i < kx.size(); ++i) {
        k.emplace_back(data, kx[i], ky, kz);
        eqForcing.emplace_back(data, k[i]);
        eqFree.emplace_back(data, k[i]);
        C.emplace_back(ZeroMatrix(4, 4));
        tMax.push_back((kxFMax - kx[i]) / ky / data.q);
    }

    double t = 0;
    spectraOut(k, C, data, t);
    addPointOfSingleSFH(k[iKFirst], C[iKFirst], data, t);

    const double tCalc  = 30;
    const double dtCalc = 5;
    for (int j = 1; j <= tCalc / dtCalc; ++j) {
        while (t < dtCalc * j) {
            #pragma omp parallel for
            for (uint i = 0; i < k.size(); ++i) {
                double tLocal   = t;
                double tMaxStep = std::min(t + dt, tMax[i]);

                bool finished = false;
                if ((k[i].x(tMaxStep) > kxFMin) && (k[i].x(t) < kxFMax)) {
                    while (finished == false) {
                        finished = eqForcing[i].make_step_forward(C[i], tLocal, tMaxStep);
                    }
                }

                finished = false;
                while (finished == false) {
                    finished = eqFree[i].make_step_forward(C[i], tLocal, t + dt);
                }
            }
            t += dt;
        }
        spectraOut(k, C, data, t);
        addPointOfSingleSFH(k[iKFirst], C[iKFirst], data, t);
    }
    return 0;
}
