// Copyright 2019 Dmitry N. Razdoburdin.

/* This file is part of DOKFUSF. DOKFUSF is a program that calculats dynamic of Keplerian flow under stochastic forcing.
DOKFUSF is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later version. DOKFUSF is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details. You should have received a copy of the GNU General Public License along with DOKFUSF; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA */

#pragma once

#include <cmath>

#include "Parameters.h"

class WaveVector {
 private:
    const double _q;

    const double _kx;
    const double _ky;
    const double _kz;

    WaveVector(double q, double kx, double ky, double kz) :
        _q(q),
        _kx(kx),
        _ky(ky),
        _kz(kz) {}

 public:
    WaveVector(const Parameters& data, double kx, double ky, double kz) :
        _q(data.q),
        _kx(kx),
        _ky(ky),
        _kz(kz) {}

    inline double x(double t) const {
        return _kx + _q * _ky * t;
    }

    inline double x() const {
        return _kx;
    }

    inline double y() const {
        return _ky;
    }

    inline double z() const {
        return _kz;
    }

    WaveVector operator() (double t) const {
        return WaveVector(_q, x(t), y(), z());
    }
};

inline double norm(const WaveVector& k) {
    return k.x() * k.x() + k.y() * k.y() + k.z() * k.z();
}

inline double abs(const WaveVector& k) {
    return sqrt(norm(k));
}

inline double norm2D(const WaveVector& k) {
    return k.x() * k.x() + k.y() * k.y();
}

inline double abs2D(const WaveVector& k) {
    return sqrt(norm2D(k));
}


