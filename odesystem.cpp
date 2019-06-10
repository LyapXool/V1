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

#include <armadillo>
#include "odesystem.hpp"
#include "instructions.hpp"

using namespace arma;
using namespace glovar;





void odesystem(const glovar::Problema problem, const bool normal, Row<double> const &x, Row<double> &f)
{
    
    switch (problem){
        case Problema::TWOORBITS:
            
            
            //Proposal
            f(0)=-1.0*x(0)*((x(0)*x(0))+(x(1)*x(1))-(1.0/4.0))*((x(0)*x(0))+(x(1)*x(1))-1.0)-x(1);
            f(1)=-1.0*x(1)*((x(0)*x(0))+(x(1)*x(1))-(1.0/4.0))*((x(0)*x(0))+(x(1)*x(1))-1.0)+x(0);
            
            break;
        case Problema::VDP:
            
            //VDP
            f(0)=x(1);
            f(1)=(1.0-(x(0)*x(0)))*x(1)-x(0);// + 1.0 / 3.0*pow(x(0), 3) - x(1);
            break;
        case Problema::HOMOCLINIC:
            
            f(0)=x(0)*(1.0-(x(0)*x(0))-(x(1)*x(1)))-x(1)*(((x(0)-1.0)*(x(0)-1.0))+(((x(0)*x(0))+x(1)-1.0)*((x(0)*x(0))+(x(1)*x(1))-1.0)));
            f(1)=x(1)*(1.0-(x(0)*x(0))-(x(1)*x(1)))+x(0)*(((x(0)-1.0)*(x(0)-1.0))+(((x(0)*x(0))+x(1)-1.0)*((x(0)*x(0))+(x(1)*x(1))-1.0)));
            break;
        case Problema::DECREASING:
            
            f(0)=1.0-x(0)*x(0);
            f(1)=-x(0)*x(1);
            break;
        case Problema::TD1:
            
            //Peter's book
            //Pag 143
            f(0) = x(0)*(pow(x(0), 2) + pow(x(1), 2) - 1.0) - x(1)*(pow(x(2), 2) + 1.0);
            f(1) = x(1)*(pow(x(0), 2) + pow(x(1), 2) - 1.0) + x(0)*(pow(x(2), 2) + 1.0);
            f(2) = 10.0 * x(2)*(pow(x(2), 2) - 1.0);
            break;
        case Problema::TD2:
            
            //On the determination of the basin of attraction of periodic orbits in three- and higher-dimensional systems
            //Journal of Mathematical Analysis and Applications
            //2009
            f(0) = -x(0)*(pow(x(0), 2) + pow(x(1), 2) - 1.0) - x(1);
            f(1) = -x(1)*(pow(x(0), 2) + pow(x(1), 2) - 1.0) + x(0);
            f(2) = -x(2);
            break;
        case Problema::TD3:
            
            //Necessary conditions for a limit cycle and its basin of attraction
            //Nonlinear Analysis 56 (2004) 643–677
            f(0) = x(0)*(-pow(x(0), 2) - pow(x(1), 2) + 1.0)*(x(0)+0.5)-x(1);
            f(1) = x(1)*(-pow(x(0), 2) - pow(x(1), 2) + 1.0)*(x(0)+0.5)+x(0);
            f(2) = 10.0*(x(2)*(x(0)+0.55)*((x(0)-0.7)*(x(0)-0.7)+0.1));
            break;
        case Problema::SIMPLE3D:
            
            //All attractor to zero
            f(0) =-x(0);
            f(1) =-x(1);
            f(2) =-x(2);
            break;
        case Problema::LORENZ:
        {
            double sx = 20.0;
            double sy = 100.0;
            double sz = 100.0;
            double si = 10.0;
            double  r = 28.0;
            double  b = 8.0 / 3.0;
            double th = -5.0 / 180.0*2.0*datum::pi;
            Mat<double> R(3,3), RI(3,3);
            Row<double> y(3);
            
            R << cos(th) <<  -sin(th) <<  0.0 << endr
            << sin(th) <<   cos(th) <<  0.0 << endr
            << 0.0     <<       0.0 <<  1.0 << endr;
            
            
            RI << cos(th) << sin(th) << 0.0 << endr
            << -sin(th)<< cos(th) << 0.0 << endr
            << 0.0     <<     0.0 << 1.0 << endr;
            
            y=x*RI;
            Row<double> fx(3);
            
            
            fx(0) = si*(-y(0) + sy / sx*y(1));
            fx(1) = r*sx / sy*y(0) - y(1) - sx*sz / sy*y(0)*y(2);
            fx(2) = -b*y(2) + sx*sy / sz*y(0)*y(1);
            f=fx*R;
            break;
        }
    }
    if(normal)
    {
        double norm2=0.0;
        double delta=10e-8;
        norm2=dot(f,f);
        f/=std::abs(sqrt(delta+norm2));
        
    }
    
    glovar::functionodecalls+=1;
}
