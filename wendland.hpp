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
 -> Authors: Carlos Argáez, Jean-Claude Berthet, Hjörtur Björnsson, Peter Giesl, Sigurdur Freyr Hafstein
 */

#ifndef wendland_hpp
#define wendland_hpp

#include <stdio.h>
#include <armadillo>

void wendlandfunction(double l,double k,double c, bool b, bool d, arma::Row<double> &wdlf);

void wendlandderivative(arma::Row<double> &wdlfinput, arma::Row<double> &wdlf1,bool negindex);

void evawdlfn(double c, double r, arma::Row<double> const &wdlfn, double &wdlfvalue,bool negindex);

#endif /* wendland_hpp */
