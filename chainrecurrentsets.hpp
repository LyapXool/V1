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

#ifndef chainrecurrentsets_hpp
#define chainrecurrentsets_hpp

#include <stdio.h>
#include <list>
#include "instructions.hpp"

void bigflowerh(const std::string gridtoeval, glovar::Problema problem, bool normal, int angles, int circles, bool spherical, double maxmaxx, double minminx, double maxmaxy, double minminy, double maxmaxz, double minminz, double alpha, double radio,arma::Mat<double> &RBFbasis,arma::Mat<double> &inputinterior, arma::Mat<double> &bigflowergrid, arma::Mat<double> &biganalgrid, std::vector<bool> &boolbigflowergrid, std::vector<bool> &boolbiganalgrid, arma::Mat<double> &cleanbigag, arma::Mat<double> &cleanbigfg, int &stride);

void qpzerpctotalh(int iteration, int angles, int circles, int &stride, const int &defcase, double critval, arma::Row<double> &qfunc, arma::Row<double> &qmatderod, arma::Row<double> &failingq, arma::Row<double> &failingqprime,  arma::Mat<double> &angrid, arma::Mat<double> &inputinterior, arma::Col<double> &newalpha, arma::Mat<double> &failinggrid, arma::Mat<double> &angridnonzero, arma::Mat<double> &flowernozero, arma::Mat<double> &collocationzero, arma::Mat<double> &collocationnozero, arma::Row<double> &alphazero, arma::Row<double> &alphanozero,std::list<int> &counter, std::list<int> &counterncol, std::vector<bool> &boolbigflowergrid);

#endif /* chainrecurrentsets_hpp */
