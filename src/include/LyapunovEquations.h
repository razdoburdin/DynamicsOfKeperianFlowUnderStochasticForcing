// Copyright 2019 Dmitry N. Razdoburdin.

/* This file is part of DOKFUSF. DOKFUSF is a program that calculats dynamic of Keplerian flow under stochastic forcing.
DOKFUSF is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later version. DOKFUSF is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details. You should have received a copy of the GNU General Public License along with DOKFUSF; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA */

#pragma once

#include <boost/numeric/ublas/matrix.hpp>

#include "Parameters.h"
#include "WaveVector.h"

namespace ublas = boost::numeric::ublas;

typedef ublas::matrix <double> Matrix;
typedef ublas::zero_matrix <double> ZeroMatrix;


inline double trace(const Matrix& C) {
    return C(0, 0) + C(1, 1) + C(2, 2) + C(3, 3);
}

inline double get_flux(const Matrix& C) {
    return C(0, 1);
}

class AbstractLyapunovEquation {
 private:
    const double _q;
    const double _invRe;
    const double _invRe_b;
    const double _Ct;

    Matrix A(double t) const;

    Matrix Adag(double t) const;

    virtual Matrix FFdag(double t) const = 0;

 protected:
    const WaveVector _k;

    static constexpr int x = 0;
    static constexpr int y = 1;
    static constexpr int z = 2;
    static constexpr int w = 3;

 public:
    double get_dt(double t) const;

    AbstractLyapunovEquation(const Parameters& data, const WaveVector& k) :
        _q(data.q),
        _invRe(data.invRe),
        _invRe_b(data.invRe_b + data.invRe / 3.0),
        _Ct(data.Ct),
        _k(k) {}

    void operator ()(const Matrix&, Matrix&, double);

    virtual void make_step_forward(Matrix&, double&) const = 0;

    inline double forsingPower(double t) const {
        return trace(FFdag(t));
    }

    virtual ~AbstractLyapunovEquation() = default;
};

class LyapunovEquationWithFlatForcing : public AbstractLyapunovEquation {
 private:
    Matrix FFdag(double) const override;

 public:
    explicit LyapunovEquationWithFlatForcing(const Parameters& data, const WaveVector& k) :
        AbstractLyapunovEquation(data, k) {}

    void make_step_forward(Matrix&, double&) const override;

    bool make_step_forward(Matrix&, double&, double) const;
};

class LyapunovEquationWith2DFlatForcing : public AbstractLyapunovEquation {
 private:
    Matrix FFdag(double) const override;

 public:
    explicit LyapunovEquationWith2DFlatForcing(const Parameters& data, const WaveVector& k) :
        AbstractLyapunovEquation(data, k) {}

    void make_step_forward(Matrix&, double&) const override;
};

class LyapunovEquationWith2DWhiteForcing : public AbstractLyapunovEquation {
 private:
    Matrix FFdag(double t) const override;

 public:
    explicit LyapunovEquationWith2DWhiteForcing(const Parameters& data, const WaveVector& k) :
        AbstractLyapunovEquation(data, k) {}

    void make_step_forward(Matrix&, double&) const override;
};

class LyapunovEquationWith3DWhiteForcing : public AbstractLyapunovEquation {
 private:
    Matrix FFdag(double t) const override;

 public:
    explicit LyapunovEquationWith3DWhiteForcing(const Parameters& data, const WaveVector& k) :
        AbstractLyapunovEquation(data, k) {}

    void make_step_forward(Matrix&, double&) const override;
};

class LyapunovEquationWith2DVorticalWhiteForcing : public AbstractLyapunovEquation {
 private:
    Matrix FFdag(double t) const override;

 public:
    explicit LyapunovEquationWith2DVorticalWhiteForcing(const Parameters& data, const WaveVector& k) :
        AbstractLyapunovEquation(data, k) {}

    void make_step_forward(Matrix&, double&) const override;
};

class LyapunovEquationWith2DSoundWhiteForcing : public AbstractLyapunovEquation {
 private:
    Matrix FFdag(double t) const override;

 public:
    explicit LyapunovEquationWith2DSoundWhiteForcing(const Parameters& data, const WaveVector& k) :
        AbstractLyapunovEquation(data, k) {}

    void make_step_forward(Matrix&, double&) const override;
};

class LyapunovEquationWithoutForcing : public AbstractLyapunovEquation {
 private:
    Matrix FFdag(double) const override;

 public:
    explicit LyapunovEquationWithoutForcing(const Parameters& data, const WaveVector& k) :
        AbstractLyapunovEquation(data, k) {}

    void make_step_forward(Matrix&, double&) const override;

    bool make_step_forward(Matrix&, double&, double) const;
};
