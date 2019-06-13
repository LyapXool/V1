/* LyapXool, version 1, is a program to compute Complete Lyapunov functions,
 -> for dynamical systems described by non linear autonomous ordinary differential equations,
 ->
 ->
 -> This program is free software; you can redistribute it and/or
 -> modify it under the terms of the GNU General Public License
 -> as published by the Free Software Foundation; either version 3
 -> of the License, or (at your option) any later version.
 ->
 -> This program is distributed in the hope that it will be useful,
 -> but WITHOUT ANY WARRANTY; without even the implied warranty of
 -> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 -> GNU General Public License for more details.
 ->
 -> You should have received a copy of the GNU General Public License
 -> along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ->
 ->
 */

#ifndef generalities_hpp
#define generalities_hpp
#include <stdio.h>
#include "instructions.hpp"

void jacobian(const glovar::Problema problem, const bool normal, arma::Row<double> const &x, arma::Mat<double> &J, const double finitedifferencetol);

void eigvalsol(const glovar::Problema problem, const bool normal, arma::Row<double> const &x, arma::cx_vec &eigval, arma::cx_mat &eigvec, const double finitedifferencetol);

void eigvecsol(const glovar::Problema problem, const bool normal, arma::Row<double> const &x, const double finitedifferencetol);

void judge(const glovar::Problema problem, bool normal, arma::Row<double> const &x, const double finitedifferencetol);

void printinformation(const bool cartesian, const std::string gridtoeval, const bool threedimensionalflower, const bool spherical, const bool normal, const int totaliterations, const double radio, const glovar::Problema problem, const double l, const double k, const double c, const double alpha, const double maxmaxx, const double minminx, const double maxmaxy, const double minminy, const double maxmaxz, const double minminz, const int circles, const int angles, const double critval);

void printhour(const int &definition);

void mcmp(long long unsigned value1, long long unsigned value2, long long unsigned &valueout);

void gcmp(long long unsigned value1, long long unsigned value2, long long unsigned &valueout);

void printall(const std::string nombre, const std::string fextension, int ordernum, arma::Mat<double> &vectoraimprimmir);

#endif /* generalities_hpp */
