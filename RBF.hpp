/* LyapXool, version 1, is a program to compute Complete Lyapunov functions,
 -> for dynamical systems described by non linear autonomous ordinary differential equations,
 ->
 -> Copyright (C) 2019 Carlos Argáez, Sigurdur Freyr Hafstein.
 ->
 -> This program is free software; you can redistribute it and/or
 -> modify it under the terms of the GNU General Public License
 -> as published by the Free Software Foundation; either version 2
 -> of the License, or (at your option) any later version.
 ->
 -> This program is distributed in the hope that it will be useful,
 -> but WITHOUT ANY WARRANTY; without even the implied warranty of
 -> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 -> GNU General Public License for more details.
 ->
 -> You should have received a copy of the GNU General Public License
 -> along with this program; if not, write carlos@hi.is or the Free Software
 -> Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 ->
 ->
 */

#ifndef RBF_hpp
#define RBF_hpp

#include <stdio.h>
#include <armadillo>
#include "instructions.hpp"

void wbase(int dimension, bool cartesian, arma::Mat<double> &wgrid);

void grid(int dimension, int maxneg, int maxpos, arma::Mat<int> &coord);

void rbfgrid(double alpha, arma::Mat<int> &gpoint, arma::Mat<double> &basis, arma::Mat<double> &outputm);

void effectivegrid(double efmaxx, double efminx, double efmaxy, double efminy, double efmaxz, double efminz,  arma::Mat<double> &gridtobeclean, arma::Mat<double> &grideffective);

void alphafunction(int dimalpha, arma::Col<double> &alp);

void alphafunctionvariable(glovar::Problema problem, bool normal, arma::Mat<double> &inputinternal, arma::Col<double> &alp);

void interpolationmatrixA(glovar::Problema problem, bool normal, double c, arma::Mat<double> &inputintec, arma::Row<double> &wedlf, arma::mat &Amat);

void choldecom(arma::mat &Amat, arma::mat &R);

void orbitalderivative(glovar::Problema problem, int iteration, bool normal, double c, arma::mat &R, arma::Col<double> &alpzro, arma::vec &betaod, arma::Mat<double> &inputinterior, arma::Row<double> &wedlf, arma::Mat<double> &domaingrid, arma::Row<double> &qmatod, arma::Row<double> &qmatderod);

#endif /* RBF_hpp */
