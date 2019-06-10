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

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <iterator>
#include <armadillo>
#include <algorithm>
#include "generalities.hpp"
#include "wendland.hpp"
#include "instructions.hpp"


using namespace std;
using namespace arma;



void wendlandfunction(double l,double k,double c, bool b, bool d, Row<double> &wdlf)
{
    wall_clock timer;
    timer.tic();
    
    long long unsigned comparing;
    long long unsigned mcm = 0;
    long long unsigned mcd = 0;
    long long unsigned N=l+1;
    long long unsigned valueout;
    Row <double>coeff((int)(N+2*k+1)), cpower((int)(N+2*k+1)), pascalcoeff((int)(N+2*k+1));
    Row <double>integral((int)(N+2*k+1)),integralbackup((int)(N+2*k+1));
    Row < long long unsigned> divisors((int)(N+2*k+1)),numerators(N+2*k+1),momentaneo((int)(N+2*k+1));
    Row <double>product((int)(N+2*k));
    wdlf.resize((int)(N+2*k+1));
    
    numerators.fill(1.0);
    
    
    for(int j=0;j<=k;++j)
    {
        product.resize((int)(N+2*j-1));
        integral.resize((int)(N+2*j));
        integralbackup.resize((int)(N+2*j));
        divisors.resize((int)(N+2*j));
        divisors.fill(1.0);
        product.zeros();
        integral.zeros();
        cpower.zeros();
        if(j==0)
        {
            if(k>0.0)
            {
                wdlf.resize(l+1);
                
                for(int i=0;i<=l;++i)
                {
                    integral.resize(l+1);
                    integralbackup.resize(l+1);
                    cpower(i)=i;
                    double x=1;
                    
                    for(int h=0;h<=i;h++)
                    {
                        coeff(h)=pow(-1,h)*x;
                        x = x * (i - h) / (h + 1);
                    }
                }
                integral=coeff;
                for(int i=0; i<l+1; ++i)
                {
                    numerators(numerators.size()-l+i-1)=coeff(i);
                }
            }else{
                for(int i=0;i<=l;++i)
                {
                    integral.resize(l+1);
                    integralbackup.resize(l+1);
                    cpower(i)=i;
                    double x=1;
                    for(int h=0;h<=i;h++)
                    {
                        coeff(h)=pow(c,cpower(h))*pow(-1,h)*x;
                        x = x * (i - h) / (h + 1);
                    }
                }
                integral=coeff;
            }
            momentaneo.fill(1.0);
        }if((j>0) && (j<k)){
            for(int i=0; i<(int)product.size()-1; ++i)
            {
                product(i+1)=coeff(i);
            }
            for(int i=1; i<(int)integral.size()-1; ++i)
            {
                integral(i+1)=product(i)/(i+1);
            }
            for(int i=(int)integral.size()-1; i>=2; i--)
            {
                divisors(i)=(i)*momentaneo(i-2);
            }
            mcm=divisors(0);
            for(int i=0; i<=(int)divisors.size()-1; ++i)
            {
                
                if((mcm == 0) || (divisors(i) == 0)){
                    break;
                }else{
                    mcmp(mcm, divisors(i),valueout);
                    mcm=(mcm * divisors(i)) / valueout; //
                }
            }
            
            divisors(0)=mcm;
            for(int i=1; i<(int)integral.size(); ++i)
            {
                integral(0)-=integral(i);
            }
            integral=-1.0*integral;
            numerators(N+2*k-l-2-j)=round(integral(0)*divisors(0));
            coeff.resize((int)integral.size());
            wdlf.resize((int)integral.size());
            momentaneo.resize((int)integral.size());
            momentaneo=divisors;
            coeff=integral;
        }if((j==k) && (k!=0)){
            for(int i=0; i<(int)product.size()-1; ++i)
            {
                product(i+1)=coeff(i);
            }
            for(int i=1; i<(int)integral.size()-1; ++i)
            {
                integral(i+1)=pow(c,i+1)*product(i)/(i+1);
            }
            for(int i=(int)integral.size()-1; i>=2; i--)
            {
                divisors(i)=(i)*momentaneo(i-2);
            }
            mcm=divisors(0);
            
            for(int i=0; i<=(int)divisors.size()-1; ++i)
            {
                if((mcm == 0) || (divisors(i) == 0)){
                    break;
                }else{
                    mcmp(mcm, divisors(i),valueout);
                    mcm=(mcm * divisors(i)) / valueout; //
                }
            }
            divisors(0)=mcm;
            for(int i=1; i<(int)integral.size(); ++i)
            {
                integral(0)-=integral(i)/pow(c,i);
            }
            integral=-1.0*integral;
            numerators(N+2*k-l-2-j)=round(integral(0)*divisors(0));
            if(numerators(N+2*k-l-2-j)==0)
            {
                numerators(N+2*k-l-2-j)=1;
            }
            for(int i=0; i<(int)integral.size()-1;++i)
            {
                if(integral(i)==0.0)
                {
                    numerators(i)=integral(i);
                }
            }
            long long unsigned mcm=divisors(0);
            
            for(int i=1; i<=(int)divisors.size()-1; ++i)
            {
                if((mcm == 0) || (divisors(i) == 0)){
                    break;
                }else{
                    mcmp(mcm, divisors(i),valueout);
                    mcm=(mcm * divisors(i)) / valueout; //
                }
            }
            
            divisors(0)=mcm;
            numerators(0)=round(divisors(0)*integral(0));
            for(int i=0; i<=(int)divisors.size()-1; ++i)
            {
                gcmp(numerators(i), divisors(i),valueout);
                mcd=valueout;
                divisors(i)=divisors(i)/mcd;
            }
            mcm=divisors(0);
            for(int i=1; i<=(int)divisors.size()-1; ++i)
            {
                if((mcm == 0) || (divisors(i) == 0)){
                    break;
                }else{
                    mcmp(mcm, divisors(i),valueout);
                    mcm=(mcm * divisors(i)) / valueout; //
                }
            }
            coeff.resize((int)integral.size());
            wdlf.resize((int)integral.size());
            momentaneo.resize((int)integral.size());
            momentaneo=divisors;
            coeff=integral;
        }
    }
    
    
    if(k==0)
    {
        wdlf=integral;
    }else{
        integralbackup=integral;
        if(b)
        {
            for(int i=0; i<(int)integral.size();++i)
            {
                wdlf(i)=round(mcm*integral(i));
            }
            Row <double> wdlfm((int)integral.size());
            wdlfm=wdlf;
            long long unsigned mcd=wdlfm(0);
            for(int i=1; i<(int)integral.size();++i)
            {
                if(wdlf(i)!=0.0)
                {
                    comparing=abs(wdlf(i));
                    gcmp(mcd, comparing,valueout);
                    mcd=valueout;
                }
            }
            if(d)
            {
                for(int i=0; i<(int)integral.size();++i)
                {
                    wdlf(i)=(wdlfm(i)/mcd);
                }
            }else{
                for(int i=0; i<(int)integral.size();++i)
                {
                    wdlf(i)=mcm*(integralbackup(i)/mcd);
                }
            }
        }else{
            wdlf=integral;
        }
    }
    
    
    glovar::outputf  << " " << endl;
    glovar::outputf  << "Computing the Wendland Function lasted: " << timer.toc() << endl;
    glovar::outputf  << " " << endl;
    
    
}


















void wendlandderivative(Row<double> &wdlfinput,Row<double> &wdlf1, bool negindex)//,Row<double> &wdlf2)
{
    negindex=false;
    int dim=(int)wdlfinput.n_cols;
    
    Row<double> temp(dim-1);
    temp.zeros();
    
    for(int i=dim-1; i>0; i--)
    {
        temp(i-1)=i*wdlfinput(i);
        
    }
    
    
    if(temp(0)!=0.0)
    {
        negindex=true;
    }
    
    if(negindex==true)
    {
        wdlf1.resize(dim-1);
        for(int i=dim-1; i>0; i--)
        {
            wdlf1(i-1)=i*wdlfinput(i);
            
        }
        
    }else{
        
        wdlf1.resize(dim-2);
        
        
        for(int i=dim-1; i-1>0; i--)
        {
            wdlf1(i-2)=i*wdlfinput(i);
            
        }
    }
    
}


void evawdlfn(double c, double r, Row<double> const &wdlfn, double &wdlfvalue, bool negindex)
{
    double battlewdlfvalue=0.0, checking;
    wdlfvalue=0.0;
    int dim = (int)wdlfn.size();
    if(negindex==false)
    {
        checking=1.0-c*r;
        if(checking>0.0)
        {
            for(int i=dim-1; i>=0; i--)
            {
                battlewdlfvalue*=r;
                battlewdlfvalue+=wdlfn(i);
            }
            wdlfvalue=battlewdlfvalue;
        }else{
            wdlfvalue=0.0;
        }
    }else{
        checking=1.0-c*r;
        if(checking>0.0)
        {
            for(int i=dim-1; i>=1; i--)
            {
                battlewdlfvalue*=r;
                battlewdlfvalue+=wdlfn(i);
            }
            battlewdlfvalue+=wdlfn(0)/r;
            wdlfvalue=battlewdlfvalue;
        }else{
            wdlfvalue=0.0;
        }
    }
}


