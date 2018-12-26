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

#include <boost/format.hpp>

#include "include/Parameters.h"
#include "include/integrator.h"
#include "include/WaveVector.h"


int main(int ac, char **av) {
    Parameters data(ac, av);
    data.output();

    const double ky    =  0.75;

    const double dk    =  0.02;
    const double kzMin =  0;
    const double kzMax =  4;
    const double kxMin = -30;
    const double kxMax =  30;

    int Nx = static_cast <int> ((kxMax - kxMin) / dk);
    int Nz = static_cast <int> ((kzMax - kzMin) / dk) + 1;

    std::stringstream mapName;
    mapName << boost::format("./map/Map(kx,kz) R = %.0le R_b = %.0le dk = %.2lf ky = %.2lf") \
        % (1.0 / data.invRe) % (1.0 / data.invRe_b) % dk % ky;
    std::ofstream fOut;
    fOut.open(mapName.str());

    omp_set_num_threads(data.Nt);
    std::vector <IntegratorOut> iOuts(Nx);
    for (int nz = 0; nz < Nz; ++nz) {
        const double kz = nz * dk + kzMin;
        #pragma omp parallel for schedule(dynamic)
        for (int nx = 0; nx < Nx; ++nx) {
            const double kx = nx * dk + kxMin;
            iOuts[nx] = integrateOverX(data, kx, kx + dk, ky, kz);
        }
        for (int nx = 0; nx < Nx; ++nx) {
            const double kx = nx * dk + kxMin;
            fOut << kx             << "\t" \
                 << kz             << "\t" \
                 << iOuts[nx].Ex   << "\t" \
                 << iOuts[nx].Ix   << "\t" \
                 << iOuts[nx].EInx << "\n";
        }
        fOut << std::endl;
    }
    fOut.close();
    return 0;
}
