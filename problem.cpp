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

#include <iostream>
#include <fstream>
#include "odesystem.hpp"
#include "instructions.hpp"
#include "wendland.hpp"
#include "RBF.hpp"
#include "chainrecurrentsets.hpp"
#include "generalities.hpp"
#include <armadillo>
#include <list>

using namespace std;
using namespace arma;
using namespace glovar;

unsigned long long int glovar::functionodecalls=0;
std::ofstream glovar::outputf;
arma::span const All=span::all;

int main()
{
   printhour(0);
    
    wall_clock timer,timer1;
    timer.tic();
    outputf.open("output.lpx", fstream::out);
    
    bool threedimensionalflower=false;
    
    int length=0;
    int stridef=0;
    
    vector<bool> boolbigflowergrid, boolbiganalgrid;
    Mat<double> criticalpoints;
    Row<double> failingqprime, failingq;
    Row<double> wdlfunction;
    Col<double> alphavector, newalpha1, alphawinter;
    Mat<double> rbfbasis;
    
    Mat<double> Cmat, bigflowergrid,biganalgrid;
    Mat<double> evalgridpoints, flowerzero, flowernozero, collocationnozero, collocationzero;
    Row<double> alphanozero,qmatod,alphazero,eigenvaluesv,qmatderod ;
    Mat<double> failinggrid, cleanbigfg, cleanbigag, griddirectional, allpos, allneg, oneone, allzero, oneposonez, onenegonez, eigenvectorsm;
    list<int> counter, counterncol;
    
    
   
    
    mat R;

    vec betaod;
   /* %%%%%%% BEGINS SECTION TO DEFINE THE CRITICAL POINTS %%%%%%%% */
    
    
    criticalpoints.resize(ncritical,variable);
    criticalpoints<<0.0 << 0.0 << endr
    <<0.0 << 0.5 << endr
    <<0.0 << 1.0 << endr;
    
    
 /*    %%%%%%% ENDS SECTION TO DEFINE THE CRITICAL POINTS %%%%%%%% */
    
    
    length=pow((abs(maxnegative)+abs(maxpositive)+1),variable);
    
    Mat<double> collocationpoints, angridnonzero, Hessian;
    
    if(variable==3)
    {
        threedimensionalflower=true;                  //is the directional a three dimensional?
    } else if(variable==2)
    {
        threedimensionalflower=false;                  //is the directional a three dimensional?
    }else{
        outputf  << "ERROR: DIMENSION. A higher dimension has not yet been implemented!!!" << endl;
        exit(9);
    }
    
    
    
    Row<double> result(variable);
    
    
    if(printing)
    {
        printinformation(cartesian, gridtoeval, threedimensionalflower, spherical, normal, totaliterations, radio, problem, l, k, c, alpha, maxmaxx, minminx, maxmaxy, minminy, maxmaxz, minminz, circles, angles, critval);
    }
    
    if(eigenvaluesjudge)
    {
        double finitedifferencetol=1e-8;
        int dim=(int)criticalpoints.n_rows;
        int dim2=(int)criticalpoints.n_cols;
        Row<double> inforvec(dim2);
        Mat<double> J(dim2,dim2);
        Mat<double> eigenvaluesm(dim2,dim2);
        cx_vec eigval;
        cx_mat eigvec;
        for(int i=0; i<dim; ++i)
        {
            inforvec=criticalpoints(i,span());
            outputf  << "The point is: " << inforvec << endl;
            jacobian(problem, normal, inforvec, J, finitedifferencetol);
            outputf  << "The Jacobian is: " << J << endl;
            eigvalsol(problem, normal, inforvec, eigval, eigvec, finitedifferencetol);
            outputf  << "The Eigenvalues are: " << eigenvaluesm << endl;
            
            judge(problem, normal, inforvec, finitedifferencetol);
        }
    }
   /*  //Second true means it is a rounded number. */
    wendlandfunction( l, k, c, true, true, wdlfunction);
    
    
    
    wbase(variable, cartesian, rbfbasis);
    
    
    {
        
        Mat<int> coordinates;
        Mat<double> rbfpoints(length,variable);
        
        grid(variable,maxnegative,maxpositive,coordinates);
        
        rbfgrid(alpha, coordinates, rbfbasis, rbfpoints);

        effectivegrid(maxmaxx,minminx,maxmaxy,minminy,maxmaxz,minminz,rbfpoints,collocationpoints);
    }
    
    

    
    

    if(constante==true)
    {
        int dimealphavector=(int)collocationpoints.n_rows;
        alphafunction(dimealphavector, alphavector);
    }else{
        alphafunctionvariable(problem, normal, collocationpoints, alphavector);
    }
    
        bigflowerh(gridtoeval,problem, normal, angles, circles, spherical, maxmaxx, minminx, maxmaxy, minminy, maxmaxz, minminz, alpha, radio, rbfbasis, collocationpoints, bigflowergrid, biganalgrid, boolbigflowergrid, boolbiganalgrid, cleanbigag, cleanbigfg, stridef);

    outputf  << "The total amount of points to evaluate your function (the length of your domain) is: " << bigflowergrid.n_rows << endl;
    {
    mat Amat;
    interpolationmatrixA(problem, normal, c, collocationpoints, wdlfunction, Amat);
    choldecom(Amat, R);
    }

    
    int i;
    for(i=0; i<=totaliterations; i++)
    {
        outputf  << "================================= Iteration no. " << i << " =================================" << endl;
        timer1.tic();
        orbitalderivative(problem, i, normal, c, R, alphavector, betaod, collocationpoints, wdlfunction, bigflowergrid, qmatod, qmatderod);

        {
            qpzerpctotalh(i, angles, circles, stridef, defcase, critval, qmatod, qmatderod, failingq, failingqprime, biganalgrid, bigflowergrid, newalpha1, failinggrid, angridnonzero, flowernozero, collocationzero, collocationnozero, alphazero, alphanozero,counter,counterncol,boolbigflowergrid);
        }

        alphavector=newalpha1;
        outputf  << "The whole proceedure for iteration no. " << i << " lasted: " << timer1.toc() << endl;
    }
    
    outputf  << " " << endl;
    outputf  << "================================= FINALIZATION =================================" << endl;
    outputf  << " " << endl;
    outputf  << "No. of times the function was called: " << functionodecalls << endl;
    outputf  << "The whole proceedure last: " << timer.toc() << endl;
    
    printhour(1);
    outputf.close();
    return 0;
    
}




