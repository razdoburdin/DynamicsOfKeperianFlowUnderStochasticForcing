// Copyright 2019 Dmitry N. Razdoburdin.

/* This file is part of DOKFUSF. DOKFUSF is a program that calculats dynamic of Keplerian flow under stochastic forcing.
DOKFUSF is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later version. DOKFUSF is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details. You should have received a copy of the GNU General Public License along with DOKFUSF; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA */

#pragma once

#include <fstream>

#include "Parameters.h"
#include "WaveVector.h"

struct IntegratorOut {
 public:
    double Ex;
    double Ix;
    double EInx;

    constexpr IntegratorOut() : Ex(0), Ix(0), EInx(0) {}

    IntegratorOut& operator += (const IntegratorOut& rhs) {
        Ex   += rhs.Ex;
        Ix   += rhs.Ix;
        EInx += rhs.EInx;

        return *this;
    }

    IntegratorOut& operator -= (const IntegratorOut& rhs) {
        Ex   -= rhs.Ex;
        Ix   -= rhs.Ix;
        EInx -= rhs.EInx;

        return *this;
    }

    IntegratorOut& operator *= (double d) {
        Ex   *= d;
        Ix   *= d;
        EInx *= d;

        return *this;
    }

    friend std::ostream& operator << (std::ostream& os, const IntegratorOut& iOut) {
        os << iOut.Ex << "\t" << iOut.Ix << "\t" << iOut.EInx;
        return os;
    }

    friend IntegratorOut operator + (IntegratorOut lhs, const IntegratorOut& rhs) {
        lhs += rhs;
        return lhs;
    }

    friend IntegratorOut operator * (IntegratorOut lhs, double d) {
        lhs *= d;
        return lhs;
    }
};

void integrationTest(const Parameters& data, const WaveVector& k, double tEnd);

IntegratorOut integrateOverX(const Parameters& data, double kxMax, double ky, double kz);

IntegratorOut integrateOverX(const Parameters& data, double kxMin, double kxMax, double ky, double kz);

void integrate(const Parameters& data, const WaveVector& kMax, double dky, double dkz);
