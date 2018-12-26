// Copyright 2019 Dmitry N. Razdoburdin.

/* This file is part of DOKFUSF. DOKFUSF is a program that calculats dynamic of Keplerian flow under stochastic forcing.
DOKFUSF is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later version. DOKFUSF is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details. You should have received a copy of the GNU General Public License along with DOKFUSF; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA */

#include "include/integrator.h"

#include <omp.h>

#include <fstream>
#include <string>
#include <vector>
#include <limits>

#include <boost/format.hpp>

#include "include/LyapunovEquations.h"
#include "include/Parameters.h"
#include "include/WaveVector.h"

inline uint32_t get_nk(const WaveVector& k, double dk) {
    return static_cast <uint32_t> ((abs(k) - std::numeric_limits <double>::epsilon()) / dk);
}

void integrationTest(const Parameters& data, const WaveVector& k, double tEnd) {
    std::ofstream fEn;
    fEn.open(data.params2Str());

    LyapunovEquationWithoutForcing eq(data, k);

    double uy  =  k.x() / sqrt(k.y() * k.y() + k.x() * k.x());
    double ux  = -k.y() / sqrt(k.y() * k.y() + k.x() * k.x());
    Matrix C = ZeroMatrix(4, 4);
    C(0, 0) = ux * ux;
    C(1, 0) = ux * uy;
    C(0, 1) = ux * uy;
    C(1, 1) = uy * uy;

    double t  = 0;
    while (t < tEnd) {
        eq.make_step_forward(C, t);
        fEn       << k.x(t) << "\t" << trace(C) << std::endl;
    }
}

IntegratorOut integrateOverX(const Parameters& data, double kxMin, double kxMax, double ky, double kz) {
    IntegratorOut iOut;

    Matrix C = ZeroMatrix(4, 4);
    double kx = kxMin;
    const WaveVector k(data, kx, ky, kz);
    LyapunovEquationWithFlatForcing eqForcing(data, k);
//    LyapunovEquationWith2DSoundWhiteForcing eqForcing(data, k);
//    LyapunovEquationWith2DVorticalWhiteForcing eqForcing(data, k);

    double t = 0;
    double dkx = 0;
    const double tMax = (kxMax - k.x()) / ky / data.q;
    bool finished = false;
    while (finished == false) {
        double E0   = trace(C);
        double I0   = get_flux(C);
        double EIn0 = eqForcing.forsingPower(t);
        double t0   = t;

        finished = eqForcing.make_step_forward(C, t, tMax);

        dkx = data.q * k.y() * (t - t0);
        iOut.Ex   += (E0 + trace(C))    * 0.5 * dkx;
        iOut.Ix   += (I0 + get_flux(C)) * 0.5 * dkx;
        iOut.EInx += (EIn0 + eqForcing.forsingPower(t)) * 0.5 * (t - t0);
    }
    iOut.Ex += trace(C)    * 0.5 * dkx;
    iOut.Ix += get_flux(C) * 0.5 * dkx;

    LyapunovEquationWithoutForcing eq(data, k);
    dkx = data.q * k.y() * eq.get_dt(t);
    iOut.Ex += trace(C)    * 0.5 * dkx;
    iOut.Ix += get_flux(C) * 0.5 * dkx;

    finished = false;
    int nSteps = 0;
    const int nStepsMin = 10;
    while (finished == false) {
        double E0 = trace(C);
        double I0 = get_flux(C);
        double t0 = t;

        eq.make_step_forward(C, t);
        dkx = data.q * k.y() * (t - t0);
        iOut.Ex += (E0 + trace(C))    * 0.5 * dkx;
        iOut.Ix += (I0 + get_flux(C)) * 0.5 * dkx;

        finished = (k(t).x() > std::abs(kxMin)) && (nSteps > nStepsMin) && (trace(C) < 0.1 * iOut.EInx);
        ++nSteps;
    }
    iOut.Ex += trace(C)    * 0.5 * dkx;
    iOut.Ix += get_flux(C) * 0.5 * dkx;

    return iOut;
}

IntegratorOut integrateOverX(const Parameters& data, double kxMax, double ky, double kz) {
    return integrateOverX(data, -kxMax, kxMax, ky, kz);
}

void integrate(const Parameters& data, const WaveVector& kMax, double dky, double dkz) {
    const double kxMax  = kMax.x();
    const double kyMax  = kMax.y();
    const double kzMax  = kMax.z();

    std::stringstream SpName;
    SpName << boost::format("Spectra(ky, kz) R = %.0le R_b = %.0le kxMax = %.0lf kyMax = %.0lf kzMax = %.0lf") \
        % (1.0 / data.invRe) % (1.0 / data.invRe_b) % kxMax % kyMax % kzMax;
    std::ofstream fOut;
    fOut.open(SpName.str());

    std::stringstream SpYName;
    SpYName << boost::format("Spectra(kz) R = %.0le R_b = %.0le kxMax = %.0lf kyMax = %.0lf kzMax = %.0lf") \
        % (1.0 / data.invRe) % (1.0 / data.invRe_b) % kxMax % kyMax % kzMax;
    std::ofstream fOutY;
    fOutY.open(SpYName.str());

    std::stringstream SpZName;
    SpZName << boost::format("Spectra(ky) R = %.0le R_b = %.0le kxMax = %.0lf kyMax = %.0lf kzMax = %.0lf") \
        % (1.0 / data.invRe) % (1.0 / data.invRe_b) % kxMax % kyMax % kzMax;
    std::ofstream fOutZ;
    fOutZ.open(SpZName.str());

    std::stringstream SpIntegratedName;
    SpIntegratedName << \
    boost::format("IntegratedSpectra R = %.0le R_b = %.0le") % (1.0 / data.invRe) % (1.0 / data.invRe_b);
    std::ofstream fOutIntegrated;
    fOutIntegrated.open(SpIntegratedName.str(), std::ios_base::app);

    omp_set_num_threads(data.Nt);

    const int Ny = static_cast <int> (kyMax / dky) - 1;
    const int Nz = static_cast <int> (kzMax / dkz);
    std::vector <std::vector <IntegratorOut>> iOuts(Ny, std::vector <IntegratorOut> (Nz));;

    const int N = Ny * Nz;
    #pragma omp parallel for schedule(dynamic)
    for (int n = 0; n < N; ++n) {
        const int ny = n % Ny;
        const int nz = n / Ny;
        const double ky = (ny + 1) * dky;
        const double kz = nz * dkz;
        iOuts[ny][nz] = integrateOverX(data, kxMax, ky, kz);
    }

    for (int ny = 0; ny < Ny; ++ny) {
        for (int nz = 0; nz < Nz; ++nz) {
            fOut << nz * dkz << "\t" << (ny + 1) * dky << "\t" << iOuts[ny][nz] << std::endl;
        }
        fOut << "\n";
    }

    std::vector <IntegratorOut> iOutsZ(Ny);
    IntegratorOut iOutZY;
    for (int ny = 0; ny < Ny; ++ny) {
        iOutsZ[ny] = (iOuts[ny][0] + iOuts[ny][Nz - 1]) * 0.5;
        for (int nz = 1; nz < Nz - 1; ++nz) {
            iOutsZ[ny] += iOuts[ny][nz];
        }
        iOutsZ[ny] *= dkz;
        iOutZY += iOutsZ[ny];
    }
    iOutZY -= (iOutsZ[0] + iOutsZ[Ny - 1]) * 0.5;
    iOutZY *= dky;

    std::vector <IntegratorOut> iOutsY(Nz);
    IntegratorOut iOutYZ;
    for (int nz = 0; nz < Nz; ++nz) {
        iOutsY[nz] = (iOuts[0][nz] + iOuts[Ny - 1][nz]) * 0.5;
        for (int ny = 1; ny < Ny - 1; ++ny) {
            iOutsY[nz] += iOuts[ny][nz];
        }
        iOutsY[nz] *= dky;
        iOutYZ += iOutsY[nz];
    }
    iOutYZ -= (iOutsY[0] + iOutsY[Nz - 1]) * 0.5;
    iOutYZ *= dkz;

    for (int ny = 0; ny < Ny; ++ny) {
        fOutZ << (ny + 1) * dky << "\t" << iOutsZ[ny] << "\n";
    }
    for (int nz = 0; nz < Nz; ++nz) {
        fOutY << nz * dkz << "\t" << iOutsY[nz] << "\n";
    }

    fOutIntegrated << data << "\t" << kxMax << "\t" << kyMax << "\t" << kzMax << "\t" << dky << "\t" << dkz << "\t" << \
        iOutYZ << "\t" << iOutZY << "\n";
}
