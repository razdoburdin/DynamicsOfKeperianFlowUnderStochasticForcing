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

    const double ky        =  1.0;
    const double kz        =  0.0;
    const double kx        = -pow(data.q * ky / data.invRe, 1.0 / 3.0);

    const double dk    =  0.005;

    std::stringstream name;
    name << boost::format("E(R) R_b = %.0le ky = %.2lf kz = %.2lf dk = %.3lf") \
                                                                        % (1.0 / data.invRe_b) % (ky) % (kz) % dk;
    std::ofstream fOut;
    fOut.open(name.str(), std::ios_base::app);

    IntegratorOut iOut = integrateOverX(data, kx, kx + dk, ky, kz);
    fOut << 1.0 / data.invRe << "\t" \
         << iOut.Ex          << "\t" \
         << iOut.Ix          << "\t" \
         << iOut.EInx        << "\n";
    fOut.close();
    return 0;
}
