# Copyright 2019 Dmitry N. Razdoburdin.

# This file is part of DOKFUSF. DOKFUSF is a program that calculats dynamic of Keplerian flow under stochastic forcing.
# DOKFUSF is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later version. DOKFUSF is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details. You should have received a copy of the GNU General Public License along with DOKFUSF; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA 

include ./configs/params.cfg
include ./configs/CPPLINT.cfg


CXXFLAGS = -O3 -march=native -ffast-math -Wall -pedantic -Wextra -Weffc++ -Wfloat-equal -Wconversion -std=c++11 -fopenmp
LDLIBS = -lboost_program_options -fopenmp

all: ./bin/IntegrationTest ./bin/Spectra[kx] ./bin/Spectra[ky,kz] ./bin/SolutionMap[kx,ky] ./bin/SolutionMap[kx,kz] ./bin/Optimal[R] ./bin/SteadyStateTransition ./bin/SolutionMap[ky,kz] ./bin/Optimal[R]

clean:
	rm -rf ./objects
	rm -rf ./bin

cheak:
	cppcheck ./src/*.cpp ./src/include/*.h
	python3 cpplint.py $(CPPLINT_OPTIONS) ./src/*.cpp ./src/include/*.h

###

steadyStateTransition: ./bin/SteadyStateTransition ./configs/params.cfg
	time ./bin/SteadyStateTransition $(KEYS)

./bin/SteadyStateTransition: ./objects/steadyStateTransition.o ./objects/Parameters.o ./objects/LyapunovEquations.o ./objects/integrator.o
	mkdir -p ./bin
	$(CXX) $(LDFLAGS) ./objects/steadyStateTransition.o ./objects/Parameters.o ./objects/LyapunovEquations.o ./objects/integrator.o -o ./bin/SteadyStateTransition $(LDLIBS)

./objects/steadyStateTransition.o: ./src/steadyStateTransition.cpp ./src/include/Parameters.h ./src/include/WaveVector.h ./src/include/LyapunovEquations.h Makefile
	mkdir -p ./objects
	$(CXX) -c $(CXXFLAGS) ./src/steadyStateTransition.cpp -o ./objects/steadyStateTransition.o

###

optimal[R]: ./bin/Optimal[R] ./configs/params.cfg
	time ./bin/Optimal[R] $(KEYS)

./bin/Optimal[R]: ./objects/optimal[R].o ./objects/Parameters.o ./objects/LyapunovEquations.o ./objects/integrator.o
	mkdir -p ./bin
	$(CXX) $(LDFLAGS) ./objects/optimal[R].o ./objects/Parameters.o ./objects/LyapunovEquations.o ./objects/integrator.o -o ./bin/Optimal[R] $(LDLIBS)

./objects/optimal[R].o: ./src/optimal[R].cpp ./src/include/Parameters.h ./src/include/WaveVector.h ./src/include/LyapunovEquations.h Makefile
	mkdir -p ./objects
	$(CXX) -c $(CXXFLAGS) ./src/optimal[R].cpp -o ./objects/optimal[R].o

###

spectra[kx]: ./bin/Spectra[kx] ./configs/params.cfg
	time ./bin/Spectra[kx] $(KEYS)

./bin/Spectra[kx]: ./objects/spectra[kx].o ./objects/Parameters.o ./objects/LyapunovEquations.o ./objects/integrator.o
	mkdir -p ./bin
	$(CXX) $(LDFLAGS) ./objects/spectra[kx].o ./objects/Parameters.o ./objects/LyapunovEquations.o ./objects/integrator.o -o ./bin/Spectra[kx] $(LDLIBS)

./objects/spectra[kx].o: ./src/spectra[kx].cpp ./src/include/Parameters.h ./src/include/WaveVector.h ./src/include/LyapunovEquations.h Makefile
	mkdir -p ./objects
	$(CXX) -c $(CXXFLAGS) ./src/spectra[kx].cpp -o ./objects/spectra[kx].o

###

integrationTest: ./bin/IntegrationTest ./configs/params.cfg
	time ./bin/IntegrationTest $(KEYS)

./bin/IntegrationTest: ./objects/integrationTest.o ./objects/Parameters.o ./objects/LyapunovEquations.o ./objects/integrator.o
	mkdir -p ./bin
	$(CXX) $(LDFLAGS) ./objects/integrationTest.o ./objects/Parameters.o ./objects/LyapunovEquations.o ./objects/integrator.o -o ./bin/IntegrationTest $(LDLIBS)

./objects/integrationTest.o: ./src/integrationTest.cpp ./src/include/Parameters.h ./src/include/integrator.h ./src/include/WaveVector.h ./src/include/LyapunovEquations.h Makefile
	mkdir -p ./objects
	$(CXX) -c $(CXXFLAGS) ./src/integrationTest.cpp -o ./objects/integrationTest.o

###

spectra[ky,kz]: ./bin/Spectra[ky,kz] ./configs/params.cfg
	time ./bin/Spectra[ky,kz] $(KEYS)

./bin/Spectra[ky,kz]: ./objects/spectra[ky,kz].o ./objects/Parameters.o ./objects/LyapunovEquations.o ./objects/integrator.o
	mkdir -p ./bin
	$(CXX) $(LDFLAGS) ./objects/spectra[ky,kz].o ./objects/Parameters.o ./objects/LyapunovEquations.o ./objects/integrator.o -o ./bin/Spectra[ky,kz] $(LDLIBS)

./objects/spectra[ky,kz].o: ./src/spectra[ky,kz].cpp ./src/include/Parameters.h ./src/include/integrator.h ./src/include/WaveVector.h ./src/include/LyapunovEquations.h Makefile
	mkdir -p ./objects
	$(CXX) -c $(CXXFLAGS) ./src/spectra[ky,kz].cpp -o ./objects/spectra[ky,kz].o

###

solutionMap[kx,ky]: ./bin/SolutionMap[kx,ky] ./configs/params.cfg
	mkdir -p ./map
	time ./bin/SolutionMap[kx,ky] $(KEYS)

./bin/SolutionMap[kx,ky]: ./objects/solutionMap[kx,ky].o ./objects/Parameters.o ./objects/LyapunovEquations.o ./objects/integrator.o
	mkdir -p ./bin
	$(CXX) $(LDFLAGS) ./objects/solutionMap[kx,ky].o ./objects/Parameters.o ./objects/LyapunovEquations.o ./objects/integrator.o -o ./bin/SolutionMap[kx,ky] $(LDLIBS)

./objects/solutionMap[kx,ky].o: ./src/solutionMap[kx,ky].cpp ./src/include/Parameters.h ./src/include/integrator.h ./src/include/WaveVector.h ./src/include/LyapunovEquations.h Makefile
	mkdir -p ./objects
	$(CXX) -c $(CXXFLAGS) ./src/solutionMap[kx,ky].cpp -o ./objects/solutionMap[kx,ky].o

###

solutionMap[kx,kz]: ./bin/SolutionMap[kx,kz] ./configs/params.cfg
	mkdir -p ./map
	time ./bin/SolutionMap[kx,kz] $(KEYS)

./bin/SolutionMap[kx,kz]: ./objects/solutionMap[kx,kz].o ./objects/Parameters.o ./objects/LyapunovEquations.o ./objects/integrator.o
	mkdir -p ./bin
	$(CXX) $(LDFLAGS) ./objects/solutionMap[kx,kz].o ./objects/Parameters.o ./objects/LyapunovEquations.o ./objects/integrator.o -o ./bin/SolutionMap[kx,kz] $(LDLIBS)

./objects/solutionMap[kx,kz].o: ./src/solutionMap[kx,kz].cpp ./src/include/Parameters.h ./src/include/integrator.h ./src/include/WaveVector.h ./src/include/LyapunovEquations.h Makefile
	mkdir -p ./objects
	$(CXX) -c $(CXXFLAGS) ./src/solutionMap[kx,kz].cpp -o ./objects/solutionMap[kx,kz].o

###

solutionMap[ky,kz]: ./bin/SolutionMap[ky,kz] ./configs/params.cfg
	mkdir -p ./map
	time ./bin/SolutionMap[ky,kz] $(KEYS)

./bin/SolutionMap[ky,kz]: ./objects/solutionMap[ky,kz].o ./objects/Parameters.o ./objects/LyapunovEquations.o ./objects/integrator.o
	mkdir -p ./bin
	$(CXX) $(LDFLAGS) ./objects/solutionMap[ky,kz].o ./objects/Parameters.o ./objects/LyapunovEquations.o ./objects/integrator.o -o ./bin/SolutionMap[ky,kz] $(LDLIBS)

./objects/solutionMap[ky,kz].o: ./src/solutionMap[ky,kz].cpp ./src/include/Parameters.h ./src/include/integrator.h ./src/include/WaveVector.h ./src/include/LyapunovEquations.h Makefile
	mkdir -p ./objects
	$(CXX) -c $(CXXFLAGS) ./src/solutionMap[ky,kz].cpp -o ./objects/solutionMap[ky,kz].o

###

solutionRelativeMap[ky,kz]: ./bin/SolutionRelativeMap[ky,kz] ./configs/params.cfg
	time ./bin/SolutionRelativeMap[ky,kz] $(KEYS)

./bin/SolutionRelativeMap[ky,kz]: ./objects/solutionRelativeMap[ky,kz].o ./objects/Parameters.o ./objects/LyapunovEquations.o ./objects/integrator.o
	mkdir -p ./map
	mkdir -p ./bin
	$(CXX) $(LDFLAGS) ./objects/solutionRelativeMap[ky,kz].o ./objects/Parameters.o ./objects/LyapunovEquations.o ./objects/integrator.o -o ./bin/SolutionRelativeMap[ky,kz] $(LDLIBS)

./objects/solutionRelativeMap[ky,kz].o: ./src/solutionRelativeMap[ky,kz].cpp ./src/include/Parameters.h ./src/include/integrator.h ./src/include/WaveVector.h ./src/include/LyapunovEquations.h Makefile
	mkdir -p ./objects
	$(CXX) -c $(CXXFLAGS) ./src/solutionRelativeMap[ky,kz].cpp -o ./objects/solutionRelativeMap[ky,kz].o

###

./objects/Parameters.o: ./src/Parameters.cpp ./src/include/Parameters.h Makefile
	$(CXX) -c $(CXXFLAGS) ./src/Parameters.cpp -o ./objects/Parameters.o

./objects/LyapunovEquations.o: ./src/LyapunovEquations.cpp ./src/include/LyapunovEquations.h ./src/include/Parameters.h ./src/include/WaveVector.h Makefile
	$(CXX) -c $(CXXFLAGS) ./src/LyapunovEquations.cpp -o ./objects/LyapunovEquations.o

./objects/integrator.o : ./src/integrator.cpp ./src/include/LyapunovEquations.h ./src/include/Parameters.h ./src/include/WaveVector.h Makefile
	$(CXX) -c $(CXXFLAGS) ./src/integrator.cpp -o ./objects/integrator.o

