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
#include "generalities.hpp"
#include "odesystem.hpp"
#include <iomanip>

using namespace arma;
using namespace glovar;
using namespace std;

arma::span const All=span::all;
void jacobian(const Problema problem, const bool normal, Row<double> const &x, Mat<double> &J, const double finitedifferencetol)
{
    wall_clock timer;
    timer.tic();
    
    

    const int N=(int)x.n_cols;
    J.resize(N,N);
    
    Row<double> xk(N), sk(N), f(N);
    
    odesystem(problem,normal,x, f);
    
    for(int i=0; i< N; ++i)
    {
        xk=x;
        xk(i)+=finitedifferencetol;
        odesystem(problem,normal,xk,sk);
        
        J(i,All)=(sk- f)/finitedifferencetol;
    }
    outputf  << "The whole proceedure to construct the Jacobian lasted: " << timer.toc() << endl;
    
    printhour(1);

}


void eigvalsol(const Problema problem, const bool normal, Row<double> const &x, cx_vec &eigval, cx_mat &eigvec, const double finitedifferencetol)
{
    wall_clock timer;
    timer.tic();
    
    

    const int N=(int)x.size();
    
    Mat<double> jacobianm(N,N);
    jacobian(problem, normal, x,jacobianm,finitedifferencetol);
    
    
    eig_gen(eigval, eigvec, jacobianm);
    outputf  << "The whole proceedure to obtain the eigen-pairs lasted: " << timer.toc() << endl;
    
    printhour(1);

}




void judge(const Problema problem, const bool normal, Row<double> const &x, const double finitedifferencetol)
{
    wall_clock timer;
    timer.tic();

    
    
    
    
    Mat<double> judgematrix((int)x.size(),(int)x.size());
    
    cx_vec eigval;
    cx_mat eigvec;
    eigvalsol(problem, normal, x, eigval, eigvec,finitedifferencetol);
    
    if(    ((real(eigval(0))!=0.0) && ((real(eigval(1))!=0.0)) && (imag(eigval(0))!=0.0) && (imag(eigval(1))!=0.0)))
    {
        if( (real(eigval(0)) < 0.0) && (real(eigval(1)) < 0.0) )
        {
            outputf  << "The critical point " << x << " is a Stable Focus  " << endl;
        }
        else if ( (real(eigval(0)) > 0.0) && (real(eigval(1)) > 0.0) )
        {
            outputf  << "The critical point " << x << " is a Unstable Focus  " << endl;
        }
    }
    
    if((imag(eigval(0))==0) && (imag(eigval(1))==0))
    {
        if(((real(eigval(0))<0) && (real(eigval(1))>0)) || ((real(eigval(0))>0) && (real(eigval(1))<0)))
        {
            outputf  << "The critical point " << x << " is a Saddle" << endl;
        }
        if( (real(eigval(0))<0) && (real(eigval(1))<0) )
        {
            outputf  << "The critical point " << x << " is a Stable Node" << endl;
        }
        if( (real(eigval(0))>0) && (real(eigval(1))>0) )
        {
            outputf  << "The critical point " << x << " is a Unstable Node" << endl;
        }
        
        
    }
    
    if( ((real(eigval(0))==0.0) && (imag(eigval(0))==0.0)) && ((imag(eigval(0))==0.0) && (imag(eigval(1))==0.0)) )
    {
        outputf  << "Linearization failed, something is wrong!! GAME OVER!! " << endl;
        
    }
    outputf  << "The whole proceedure to Judge the critical point lasted: " << timer.toc() << endl;
    
    printhour(1);

    
}


void printinformation(const bool cartesian, const string gridtoeval, const bool threedimensionalflower, const bool spherical, const bool normal, const int totaliterations, const double radio, const Problema problem, const double l, const double k, const double c, const double alpha, const double maxmaxx, const double minminx, const double maxmaxy, const double minminy, const double maxmaxz, const double minminz, const int circles, const int angles, const double critval)
{
    ofstream datos;
    
    datos.open ("datos.lpx", fstream::app);
    
    if(cartesian==false){
        datos << "Hexagonal grid" << endl;
    }else{
        datos << "Canonical grid" << endl;
    }
    
    if(gridtoeval=="directional")
    {
        {
            datos << "It will use directional flower" << endl;
        }
    }
    if(gridtoeval=="circular"){
        if(threedimensionalflower)
        {
            datos << "It will use the directional 3Dflower ";
        }else{
            datos << "It will use the directional flower ";
        }
    }
    if(spherical)
    {
        datos << "It will use the spherical flower ";
    }else {
            datos << "It will use the flower along the axis with RBF set of vector basis ";
    }
    if(normal)
    {
        datos << "under the Delta method approach ";
    }else{
        datos << "under the unnormalized method approach ";
    }
    if(totaliterations>0)
    {
        datos << "for a total of " << totaliterations << " iterations." << endl;
    }else{
        datos << "for only the computation of the Lyapunov function." << endl;
    }
    
    datos << "The general settings for this problem are: " << endl;
    datos << "Problem: " << probnames[problem] << " | " <<  " Wendland Function Parameters:"  << " l=" << l << ", k=" << k << ", c=" << c << " || " << endl;
    datos << "alpha: " << alpha << " | " << endl;
    datos << "x-maxmax: " << maxmaxx << " | " << " x-minmin " << minminx  << " || " <<endl;
    datos << "y-maxmax: " << maxmaxy << " | " << " y-minmin " << minminy  << " || " << endl;
    datos << "z-maxmax: " << maxmaxz << " | " << " z-minmin " << minminz  << " || " << endl;
    datos << "radio/interval: " << radio << " | "<< " circles: " << circles << " || "  << endl;
    datos << "angles: " << angles << " | " << " critval: " << critval << " || " << endl;
    datos << " " << endl;
    if(cartesian)
    {
        datos << "THE COLLOCATION GRID WAS SET TO BE CARTESIAN" << endl;
    }
    datos.close();
}







void printall(const string nombre, const string fextension, int ordernum, Mat<double> &vectoraimprimmir)
{
    wall_clock timer;
    timer.tic();
    
    int dim=(int)vectoraimprimmir.n_rows;
    int dim2=(int)vectoraimprimmir.n_cols;
    
    if((dim!=0)&&(dim2!=0))
    {
    
    

    outputf << "Printing results..." << endl;

    ostringstream fileName;
    
    fileName<<"s"<<nombre<<"."<<fextension;
    ofstream valuesdocument(fileName.str(), fstream::out | fstream::app);
    
    valuesdocument << nombre << ordernum << "=[";
    if(vectoraimprimmir.n_cols > vectoraimprimmir.n_rows)
    {
        for(int p=0; p<dim2; ++p)
        {
            valuesdocument << std::fixed << std::setprecision(18) << vectoraimprimmir(0,p) << " " ;
        }
    }
    else
    {
        for(int p=0; p<dim; ++p)
        {
            valuesdocument << std::fixed << std::setprecision(18) << vectoraimprimmir(p,0) << " " ;
        }
    }
    
    valuesdocument << "];" << endl;
    valuesdocument.close();
    outputf  << "The whole proceedure to print results lasted: " << timer.toc() << endl;
    
    printhour(1);
    }else{
        outputf << "WARNING: " << nombre << " does not contain values and its dimension is 0 " << endl;
    }

}


void printhour(const int &definition)
{
    
    
    
    time_t tempus;
    struct tm * infotiempo;
    char datos[100];
    
    time (&tempus);
    infotiempo = localtime(&tempus);
    if(definition==0)
    {
        strftime(datos,sizeof(datos),"Computation started on %d-%m-%Y at %H:%M:%S",infotiempo);
        
    }if(definition==1)
    {
        strftime(datos,sizeof(datos),"Computation finished on %d-%m-%Y at %H:%M:%S",infotiempo);
    }
    string str(datos);
    
    outputf  << str << endl;
    
}

void mcmp(long long unsigned value1, long long unsigned value2, long long unsigned &valueout)
{
    while (value2 > 0) {
        long long unsigned r = value1 % value2;
        value1 = value2;
        value2 = r;
    }
    valueout=value1;
}

void gcmp(long long unsigned value1, long long unsigned value2, long long unsigned &valueout)
{
    while (value2 > 0) {
        long long unsigned r = value1 % value2;
        value1 = value2;
        value2 = r;
    }
    valueout=value1;
}


