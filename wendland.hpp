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
#ifndef wendland_hpp
#define wendland_hpp

#include <stdio.h>
#include <armadillo>

void wendlandfunction(double l,double k,double c, bool b, bool d, arma::Row<double> &wdlf);

void wendlandderivative(arma::Row<double> &wdlfinput, arma::Row<double> &wdlf1,bool negindex);

void evawdlfn(double c, double r, arma::Row<double> const &wdlfn, double &wdlfvalue,bool negindex);

#endif /* wendland_hpp */
