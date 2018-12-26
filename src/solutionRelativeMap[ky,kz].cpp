// Copyright 2019 Dmitry N. Razdoburdin.

/* This file is part of DOKFUSF. DOKFUSF is a program that calculats dynamic of Keplerian flow under stochastic forcing.
DOKFUSF is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later version. DOKFUSF is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details. You should have received a copy of the GNU General Public License along with DOKFUSF; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA */

#include <omp.h>

#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include <boost/format.hpp>

#include "include/Parameters.h"
#include "include/integrator.h"
#include "include/WaveVector.h"


int main(int ac, char **av) {
    Parameters data(ac, av);
    data.output();

    const double dk    =  0.02;

    const double kyMin =  dk;
    const double kyMax =  4;

    const double kzMin =  0;
    const double kzMax =  4;


    int Ny = static_cast <int> ((kyMax - kyMin) / dk) + 1;
    int Nz = static_cast <int> ((kzMax - kzMin) / dk) + 1;

    std::stringstream mapName;
    mapName << boost::format("./map/RelMap(ky,kz) R = %.0le R_b = %.0le dk = %.3lf kx = optimal") \
        % (1.0 / data.invRe) % (1.0 / data.invRe_b) % dk;
    std::ofstream fOut;
    fOut.open(mapName.str());

    omp_set_num_threads(data.Nt);
    std::vector <IntegratorOut> iOuts(Ny);
    std::vector <IntegratorOut> iOutsFlat(Ny);
    for (int nz = 0; nz < Nz; ++nz) {
        const double kz = nz * dk + kzMin;
        #pragma omp parallel for schedule(dynamic)
        for (int ny = 0; ny < Ny; ++ny) {
            const double ky = ny * dk + kyMin;
            const double kx = -pow(data.q / data.invRe * ky, 1.0 / 3.0);
            iOuts[ny] = integrateOverX(data, kx, kx + dk, ky, kz);
        }

        if (nz == 0) {
            iOutsFlat = iOuts;
        }

        for (int ny = 0; ny < Ny; ++ny) {
            const double ky = ny * dk + kyMin;
            fOut << ky                                  << "\t" \
                 << kz                                  << "\t" \
                 << iOuts[ny].Ex                        << "\t" \
                 << iOuts[ny].Ix                        << "\t" \
                 << iOuts[ny].EInx                      << "\t" \
                 << iOuts[ny].Ex   / iOutsFlat[ny].Ex   << "\t" \
                 << iOuts[ny].Ix   / iOutsFlat[ny].Ix   << "\t" \
                 << iOuts[ny].EInx / iOutsFlat[ny].EInx << "\n";
        }
        fOut << std::endl;
    }
    fOut.close();
    return 0;
}
