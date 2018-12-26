// Copyright 2019 Dmitry N. Razdoburdin.

/* This file is part of DOKFUSF. DOKFUSF is a program that calculats dynamic of Keplerian flow under stochastic forcing.
DOKFUSF is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later version. DOKFUSF is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details. You should have received a copy of the GNU General Public License along with DOKFUSF; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA */

#include "include/LyapunovEquations.h"

#include <algorithm>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/odeint/integrate/integrate.hpp>
#include <boost/numeric/odeint.hpp>

namespace ode = boost::numeric::odeint;

Matrix AbstractLyapunovEquation::A(double t) const {
    Matrix M = ZeroMatrix(4, 4);

//             non-viscous  viscosity                bulk viscosity
    M(x, x) =              -norm(_k(t)) * _invRe -_k.x(t) * _k.x(t) * _invRe_b;
    M(x, y) =  2                                 -_k.y()  * _k.x(t) * _invRe_b;
    M(x, z) =                                    -_k.z()  * _k.x(t) * _invRe_b;
    M(x, w) =  _k.x(t);

    M(y, x) = -(2 - _q)                          -_k.x(t) * _k.y()  * _invRe_b;
    M(y, y) =              -norm(_k(t)) * _invRe -_k.y()  * _k.y()  * _invRe_b;
    M(y, z) =                                    -_k.z()  * _k.y()  * _invRe_b;
    M(y, w) =  _k.y();

    M(z, x) =                                    -_k.x(t) * _k.z()  * _invRe_b;
    M(z, y) =                                    -_k.y()  * _k.z()  * _invRe_b;
    M(z, z) =              -norm(_k(t)) * _invRe -_k.z()  * _k.z()  * _invRe_b;
    M(z, w) =  _k.z();

    M(w, x) = -_k.x(t);
    M(w, y) = -_k.y();
    M(w, z) = -_k.z();
    M(w, w) =  0;
    return M;
}

Matrix AbstractLyapunovEquation::Adag(double t) const {
    Matrix M = ZeroMatrix(4, 4);

//             non-viscous  viscosity               bulk viscosity
    M(x, x) =              -norm(_k(t)) * _invRe -_k.x(t) * _k.x(t) * _invRe_b;
    M(x, y) = -(2 - _q)                          -_k.y()  * _k.x(t) * _invRe_b;
    M(x, z) =                                    -_k.z()  * _k.x(t) * _invRe_b;
    M(x, w) =  -_k.x(t);

    M(y, x) =  2                                 -_k.x(t) * _k.y()  * _invRe_b;
    M(y, y) =              -norm(_k(t)) * _invRe -_k.y()  * _k.y()  * _invRe_b;
    M(y, z) =                                    -_k.z()  * _k.y()  * _invRe_b;
    M(y, w) =  -_k.y();

    M(z, x) =                                    -_k.x(t) * _k.z()  * _invRe_b;
    M(z, y) =                                    -_k.y()  * _k.z()  * _invRe_b;
    M(z, z) =              -norm(_k(t)) * _invRe -_k.z()  * _k.z()  * _invRe_b;
    M(z, w) =  -_k.z();

    M(w, x) =  _k.x(t);
    M(w, y) =  _k.y();
    M(w, z) =  _k.z();
    M(w, w) =  0;
    return M;
}

double AbstractLyapunovEquation::get_dt(double t) const {
    double stepSound = _Ct / abs(_k(t));
    if (_invRe > 0) {
        double stepRe    = _Ct / (_invRe   * norm(_k(t)));
        double stepRe_b  = _Ct / (_invRe_b * norm(_k(t)));

        double step = std::min(stepSound, stepRe);
        return std::min(step, stepRe_b);
    } else {
        return stepSound;
    }
}

void AbstractLyapunovEquation::operator ()(const Matrix& C, Matrix& dCdt, double t) {
    dCdt = FFdag(t) + ublas::prod(A(t), C) + ublas::prod(C, Adag(t));
}

Matrix LyapunovEquationWithFlatForcing::FFdag(double) const {
    Matrix M = ZeroMatrix(4, 4);
    constexpr double therd = 1.0 / 3.0;
    M(x, x) = therd;
    M(y, y) = therd;
    M(z, z) = therd;
    return M;
}

void LyapunovEquationWithFlatForcing::make_step_forward(Matrix &C, double& t) const {
    double dt = get_dt(t);
    ode::integrate(*this, C, t, t + dt, dt);
    t += dt;
}

bool LyapunovEquationWithFlatForcing::make_step_forward(Matrix &C, double& t, double tMax) const {
    double dt = get_dt(t);
    if (t + dt < tMax) {
        ode::integrate(*this, C, t, t + dt, dt);
        t += dt;
        return false;
    } else {
        ode::integrate(*this, C, t, tMax, dt);
        t = tMax;
        return true;
    }
}

Matrix LyapunovEquationWith2DFlatForcing::FFdag(double) const {
    Matrix M = ZeroMatrix(4, 4);
    M(x, x) = 0.5;
    M(y, y) = 0.5;
    return M;
}

void LyapunovEquationWith2DFlatForcing::make_step_forward(Matrix &C, double& t) const {
    double dt = get_dt(t);
    ode::integrate(*this, C, t, t + dt, dt);
    t += dt;
}

Matrix LyapunovEquationWith2DWhiteForcing::FFdag(double t) const {
    Matrix M = ZeroMatrix(4, 4);
    M(x, x) = 1 / abs2D(_k(t));
    M(y, y) = 1 / abs2D(_k(t));
    M(z, z) = 1 / abs2D(_k(t));
    return M;
}

void LyapunovEquationWith2DWhiteForcing::make_step_forward(Matrix &C, double& t) const {
    double dt = get_dt(t);
    ode::integrate(*this, C, t, t + dt, dt);
    t += dt;
}

Matrix LyapunovEquationWith3DWhiteForcing::FFdag(double t) const {
    Matrix M = ZeroMatrix(4, 4);
    M(x, x) = 1 / norm(_k(t));
    M(y, y) = 1 / norm(_k(t));
    M(z, z) = 1 / norm(_k(t));
    return M;
}

void LyapunovEquationWith3DWhiteForcing::make_step_forward(Matrix &C, double& t) const {
    double dt = get_dt(t);
    ode::integrate(*this, C, t, t + dt, dt);
    t += dt;
}

Matrix LyapunovEquationWith2DVorticalWhiteForcing::FFdag(double t) const {
    Matrix M = ZeroMatrix(4, 4);
    M(x, x) =  _k.y()  * _k.y()  / norm(_k(t));
    M(x, y) = -_k.x(t) * _k.y()  / norm(_k(t));
    M(y, x) = -_k.x(t) * _k.y()  / norm(_k(t));
    M(y, y) =  _k.x(t) * _k.x(t) / norm(_k(t));
    return M;
}

void LyapunovEquationWith2DVorticalWhiteForcing::make_step_forward(Matrix &C, double& t) const {
    double dt = get_dt(t);
    ode::integrate(*this, C, t, t + dt, dt);
    t += dt;
}

Matrix LyapunovEquationWith2DSoundWhiteForcing::FFdag(double t) const {
    Matrix M = ZeroMatrix(4, 4);
    M(x, x) = _k.x(t) * _k.x(t) / norm(_k(t));
    M(x, y) = _k.x(t) * _k.y()  / norm(_k(t));
    M(y, x) = _k.x(t) * _k.y()  / norm(_k(t));
    M(y, y) = _k.y()  * _k.y()  / norm(_k(t));
    return M;
}

void LyapunovEquationWith2DSoundWhiteForcing::make_step_forward(Matrix &C, double& t) const {
    double dt = get_dt(t);
    ode::integrate(*this, C, t, t + dt, dt);
    t += dt;
}


Matrix LyapunovEquationWithoutForcing::FFdag(double) const {
    return ZeroMatrix(4, 4);
}

void LyapunovEquationWithoutForcing::make_step_forward(Matrix &C, double& t) const {
    double dt = get_dt(t);
    ode::integrate(*this, C, t, t + dt, dt);
    t += dt;
}

bool LyapunovEquationWithoutForcing::make_step_forward(Matrix &C, double& t, double tMax) const {
    double dt = get_dt(t);
    if (t + dt < tMax) {
        ode::integrate(*this, C, t, t + dt, dt);
        t += dt;
        return false;
    } else {
        ode::integrate(*this, C, t, tMax, dt);
        t = tMax;
        return true;
    }
}
