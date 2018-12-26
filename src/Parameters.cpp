// Copyright 2019 Dmitry N. Razdoburdin.

/* This file is part of DOKFUSF. DOKFUSF is a program that calculats dynamic of Keplerian flow under stochastic forcing.
DOKFUSF is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later version. DOKFUSF is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details. You should have received a copy of the GNU General Public License along with DOKFUSF; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA */

#include "include/Parameters.h"

#include <iostream>
#include <array>
#include <string>
#include <sstream>
#include <fstream>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

Parameters::ParamsArray Parameters::InitParams(int ac, char* av[]) {
    Parameters::ParamsArray pA;
    std::string Re;
    std::string Re_b;

    po::options_description data("Allowed options");
    data.add_options()
     ("help", "produce help message")
     ("q",        po::value <double> (&pA[qPosition])        -> default_value(1.5),     "Shear rate")
     ("R",        po::value <std::string> (&Re)              -> default_value("inf"),   "Reynolds number")
     ("R_b",      po::value <std::string> (&Re_b)            -> default_value("inf"),   "Second Reynolds number")
     ("Ct",       po::value <double> (&pA[CtPosition])       -> default_value(0.1),     "Courant constant")
     ("Nt",       po::value <double> (&pA[NtPosition])       -> default_value(1),       "Number of threads");

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, data), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << data << std::endl;
    }

    if (Re.compare("inf") == 0) {
//        fprintf(stderr, "%s\n%s\n", "R is set to be infinite. In that case saturation is impossible.", "Calculations will not finished any when! Program is stopped");
//        std::exit(EXIT_FAILURE);
        pA.at(invRePosition) = 0;
    } else {
        std::stringstream ss;
        double x;
        ss << Re;
        ss >> x;
        if (x < 1e-3) {
            fprintf(stdout, "%s\n", "R cannot be less than 1e-3. Re=inf is used by default");
            pA.at(invRePosition) = 0;
        } else {
            pA.at(invRePosition) = 1.0 / x;
        }
    }

    if (Re_b.compare("inf") == 0) {
        pA.at(invRe_bPosition) = 0;
    } else {
        std::stringstream ss;
        double x;
        ss << Re_b;
        ss >> x;
        if (x < 1e-3) {
            std::cout << "R_b cannot be less than 1e-3. R_b=inf is used by default" << std::endl;
            pA.at(invRe_bPosition) = 0;
        } else {
            pA.at(invRe_bPosition) = 1.0 / x;
        }
    }

    return pA;
}

Parameters::Parameters(int ac, char** av) :
    _pA(InitParams(ac, av)),
    q(_pA[qPosition]),
    invRe(_pA.at(invRePosition)),
    invRe_b(_pA.at(invRe_bPosition)),
    Ct(_pA.at(CtPosition)),
    Nt(static_cast <int> (_pA.at(NtPosition))) {}

std::string Parameters::params2Str() const {
    std::stringstream ss;
    ss << "q        = " << q             << " ";
    ss << "R        = " << 1.0 / invRe   << " ";
    ss << "R_b      = " << 1.0 / invRe_b << " ";
    ss << "Ct       = " << Ct;
    return ss.str();
}

std::ofstream& operator << (std::ofstream& os, const Parameters& data) {
    os << data.q << "\t" << 1.0 / data.invRe << "\t" << 1.0 / data.invRe_b;
    return os;
}

void Parameters::output() const {
    fprintf(stdout, "%s\n", params2Str().c_str());
}
