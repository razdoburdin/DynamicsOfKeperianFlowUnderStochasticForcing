// Copyright 2019 Dmitry N. Razdoburdin.

/* This file is part of DOKFUSF. DOKFUSF is a program that calculats dynamic of Keplerian flow under stochastic forcing.
DOKFUSF is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later version. DOKFUSF is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details. You should have received a copy of the GNU General Public License along with DOKFUSF; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA */

#include "include/Parameters.h"
#include "include/WaveVector.h"
#include "include/integrator.h"


int main(int ac, char **av) {
    Parameters data(ac, av);
    data.output();

    integrationTest(data, WaveVector(data, -20, 1, 0), 30);

    return 0;
}
