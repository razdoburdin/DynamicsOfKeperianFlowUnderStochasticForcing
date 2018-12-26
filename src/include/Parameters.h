// Copyright 2019 Dmitry N. Razdoburdin.

/* This file is part of DOKFUSF. DOKFUSF is a program that calculats dynamic of Keplerian flow under stochastic forcing.
DOKFUSF is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later version. DOKFUSF is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details. You should have received a copy of the GNU General Public License along with DOKFUSF; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA */

#pragma once

#include <array>
#include <string>
#include <cmath>
#include <fstream>

class Parameters {
 private:
    static constexpr int NParams = 5;

    static constexpr const double PI = std::atan(1.0) * 4;

    typedef std::array <double, NParams> ParamsArray;

    ParamsArray _pA;
    ParamsArray InitParams(int, char**);

    static constexpr int qPosition        = 0;
    static constexpr int invRePosition    = 1;
    static constexpr int invRe_bPosition  = 2;
    static constexpr int CtPosition       = 3;
    static constexpr int NtPosition       = 4;

 public:
    const double q;
    const double invRe;
    const double invRe_b;
    const double Ct;
    const int Nt;

    Parameters(int, char**);

    std::string params2Str() const;

    void output() const;
};

std::ofstream& operator << (std::ofstream& os, const Parameters& data);
