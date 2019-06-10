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

#include <math.h>
#include <iostream>
#include <sstream>
#include <armadillo>
#include <fstream>
#include <list>
#include <string>
#include "RBF.hpp"
#include "instructions.hpp"
#include "odesystem.hpp"
#include "wendland.hpp"
#include <numeric>
#include "generalities.hpp"
#if defined(_OPENMP)
#include <omp.h>
#endif
using namespace arma;
using namespace std;
using namespace glovar;

arma::span const All=span::all;

void wbase(int dimension, bool cartesian, Mat<double> &wgrid)
{
    wall_clock timer;
    timer.tic();
    if(cartesian)
    {
        wgrid.resize(dimension,dimension);
        wgrid.eye(dimension,dimension);
    }else{
        
        Row<double> ek(dimension);
        wgrid.resize(dimension,dimension);
        wgrid.zeros();
        ek.zeros();
        for(int k=1; k<=dimension; ++k)
        {
            
            ek(k-1)=sqrt(1.0/(2*k*(k+1)));
            wgrid(k-1,span::all)=ek;
            
            wgrid(k-1,All)(k-1)=(k+1)*ek(k-1);
        }
        
    }
    
    if(printing)
    {
        if(dimension==2)
        {
            Row<double> x=wgrid(0,All);
            Row<double> y=wgrid(1,All);
            printall("xrbfbasis",fextension,0,x);
            printall("yrbfbasis",fextension,0,y);
        }
        if(dimension==3)
        {
            Row<double> x=wgrid(0,All);
            Row<double> y=wgrid(1,All);
            Row<double> z=wgrid(2,All);
            printall("xrbfbasis",fextension,0,x);
            printall("yrbfbasis",fextension,0,y);
            printall("zrbfbasis",fextension,0,z);
            
        }
    }
    
    
    outputf  << " " << endl;
    outputf  << "Computing the base lasted: " << timer.toc() << endl;
    
    printhour(1);
}


void grid(int dimension, int maxneg, int maxpos, Mat<int> &coord)
{
    wall_clock timer;
    timer.tic();

    int i=0, j=0, k=0;
    
    int elements;
    int mn=1;
    coord.zeros();
    elements=abs(maxneg)+abs(maxpos)+1;
    vector<int> m(dimension);
    vector<vector<int>> v(dimension);
    
    {
        for(i=0; i<dimension;++i)
        {
            v[i].resize(elements);
            for(j=0;j<elements;++j)
            {
                v[i][j]=maxneg+j;
            }
        }
        
        for(i=0; i<dimension; ++i)
        {
            m[i]=(int)v[i].size();
            mn*=m[i];
        }
        coord.resize((int)mn,(int)dimension);
        
        for(i=0; i<mn;++i){
            k=i;
            for(j=dimension-1;j>=0;--j){
                coord(i,j)=v[j][k%m[j]];
                k/=m[j];
            }
        }
    }
    
    
    outputf  << " " << endl;
    outputf  << "Computing the grid lasted: " << timer.toc() << endl;
    
    printhour(1);
}


void rbfgrid(double alpha, Mat<int> &gpoint, Mat<double> &basis, Mat<double> &outputm)
{
    
    wall_clock timer;
    timer.tic();
#if defined(_OPENMP)
    omp_set_dynamic(0);
    omp_set_num_threads(OMP_NUM_THREADS);
#endif
    int i=0, j=0, tid=0;

    
    

    
    outputm.zeros();//=0.0;
    int dimgpoint1=(int)gpoint.n_rows;
    int dimgpoint2=(int)gpoint.n_cols;
    
    int chunk = int(floor(dimgpoint1/OMP_NUM_THREADS));
#pragma omp parallel num_threads(OMP_NUM_THREADS) shared(outputm,chunk) private(tid,i,j)
    {
#if defined(_OPENMP)
        tid=omp_get_thread_num();
        if(tid==0)
        {
            outputf <<"Starting parallel region\n" << endl;
        }
#endif
#pragma omp barrier
#pragma omp for schedule (dynamic, chunk) nowait
        for(i=0; i<dimgpoint1; ++i)
        {
            for(j=0; j<dimgpoint2; ++j)
            {
                outputm(i,All)+=alpha*gpoint(i,j)*basis(j,All);
            }
            if(dimgpoint2==2)
            {
                outputm(i,All)+=0.5*alpha*(basis(0,All)+basis(1,All));
            }
            if(dimgpoint2==3)
            {
                outputm(i,All)+=0.5*alpha*(basis(0,All)+basis(1,All)+basis(2,All));
            }
        }
        
    }
    bool printing2=false;
    if(printing2)
    {
        if(dimgpoint2==2)
        {
            Col<double> x=outputm(All,0);
            Col<double> y=outputm(All,1);
            printall("xrbfgrid",fextension,0,x);
            printall("yrbfgrid",fextension,0,y);
        }
        if(dimgpoint2==3)
        {
            Col<double> x=outputm(All,0);
            Col<double> y=outputm(All,1);
            Col<double> z=outputm(All,2);
            printall("xrbfgrid",fextension,0,x);
            printall("yrbfgrid",fextension,0,y);
            printall("yrbfgrid",fextension,0,z);
            
        }
    }
    outputf  << " " << endl;
    outputf  << "Computing the RBF grid lasted: " << timer.toc() << endl;
    
    printhour(1);
}


void effectivegrid(double efmaxx, double efminx, double efmaxy, double efminy, double efmaxz, double efminz, Mat<double> &gridtobeclean, Mat<double> &grideffective)
{
    wall_clock timer;
    timer.tic();
    int width=(int)gridtobeclean.n_cols;
    
    
    if(width==2)
    {
        if(efmaxx<=efminx)
        {
            outputf  << "ERROR: Maximum should be larger than Minimum" << endl;
            exit(9);
        }
        if(efmaxy<=efminy)
        {
            outputf  << "ERROR: Maximum should be ddlarger than Minimum" << endl;
            exit(9);
        }
    }
    if(width==3)
    {
        if(efmaxx<=efminx)
        {
            outputf  << "ERROR: Maximum should be larger than Minimum" << endl;
            exit(9);
        }
        if(efmaxy<=efminy)
        {
            outputf  << "ERROR: Maximum should be larger than Minimum" << endl;
            exit(9);
        }
        if((efmaxz==0.0)&&(efminz==0.0))
        {
            efmaxz=efmaxx;
            efminz=-efmaxz;
        }
        if(efmaxz<=efminz)
        {
            outputf  << "ERROR: Maximum should be larger than Minimum" << endl;
            exit(9);
        }
    }
    
    int dim1=(int)gridtobeclean.n_rows;//longitud
    int dim2=(int)gridtobeclean.n_cols;//anchura
    
    list<int> counter;
    if(dim2==2)
    {
        
        for(int i=0; i<dim1; ++i)
        {
            if((gridtobeclean(i,0)<=efmaxx && gridtobeclean(i,0)>=efminx) && (gridtobeclean(i,1)<=efmaxy && gridtobeclean(i,1)>=efminy))
            {
                counter.push_back(i);
            }
        }
    }
    
    if(dim2==3)
    {
        for(int i=0; i<dim1; ++i)
        {
            if((gridtobeclean(i,0)<=efmaxx && gridtobeclean(i,0)>=efminx) && (gridtobeclean(i,1)<=efmaxy && gridtobeclean(i,1)>=efminy) && (gridtobeclean(i,2)<=efmaxz && gridtobeclean(i,2)>=efminz))
            {
                counter.push_back(i);
            }
        }
    }
    
    int fin=(int)counter.size();
    grideffective.resize(fin,dim2);
    int n=0;
    for(list<int>::iterator i=counter.begin(); i!=counter.end(); ++i)
    {
        grideffective(n,All)=gridtobeclean(*i,All);
        n++;
    }
    
    
    
    counter.clear();
    
    
    if(glovar::printing)
    {
        int iteration=0;
        if(dim2==2)
        {
            Col<double> x(fin);
            Col<double> y(fin);
            x=grideffective(All,0);
            y=grideffective(All,1);
            printall("xcollocationpoints", fextension, iteration, x);
            printall("ycollocationpoints", fextension, iteration, y);
        }
        if(dim2==3)
        {
            Col<double> x(fin);
            Col<double> y(fin);
            Col<double> z(fin);
            x=grideffective(All,0);
            y=grideffective(All,1);
            z=grideffective(All,2);
            printall("xcollocationpoints", fextension, iteration, x);
            printall("ycollocationpoints", fextension, iteration, y);
            printall("zcollocationpoints", fextension, iteration, z);
        }
    }
    
    outputf  << " " << endl;
    outputf  << "Constraining the RBF grid to the boundaries lasted: " << timer.toc() << endl;
    
    printhour(1);
}

void alphafunction(int dimalpha, Col<double> &alp)
{
    wall_clock timer;
    timer.tic();
    
    
    alp.resize(dimalpha);
    alp.fill(-1.0);
    outputf  << " " << endl;
    outputf  << "The whole proceedure for the alpha vector lasted: " << timer.toc() << endl;
    
    printhour(1);
}


void alphafunctionvariable(glovar::Problema problem, bool normal, Mat<double> &inputinternal, Col<double> &alp)
{
    wall_clock timer;
    timer.tic();
    
    
    int dim=(int)inputinternal.n_cols;
    int dim2=(int)inputinternal.n_rows;
    alp.resize(dim2);
    Row<double> result(dim),saving(dim);
    double norm1,norm2;
    
    for(int i=0; i<dim2; ++i)
    {
        saving=inputinternal(i,All);
        norm1=dot(saving,saving);
        odesystem(problem,normal,saving, result);
        norm2=dot(result,result);
        alp(i)=-1.0*norm1*(1.0+norm2);
    }
    outputf  << " " << endl;
    outputf  << "The whole proceedure for the alpha vector variable lasted: " << timer.toc() << endl;
    
    printhour(1);
}



void interpolationmatrixA(glovar::Problema problem, bool normal, double c, Mat<double> &inputintec, Row<double> &wedlf, mat &Amat)
{
    wall_clock timer;
#if defined(_OPENMP)
    omp_set_dynamic(0);
    omp_set_num_threads(OMP_NUM_THREADS);
#endif
    int j=0, k=0, tid=0;
    
    
    if(normal)
    {
        outputf  << "Computing Interpolation Matrix with almost normalized function" << endl;
    }else{
        outputf  << "Computing Interpolation Matrix with no normalized function" << endl;
    }
    int dimwf=(int)wedlf.size();
    
    Row<double> wdlf1(dimwf-2),wdlf2(dimwf-4);
    
    double atzero;
    bool negindex=false;
    wendlandderivative(wedlf,wdlf1,negindex);
    
    
    wendlandderivative(wdlf1,wdlf2,negindex);
    
    
    int dimA=(int)inputintec.n_rows;
    
    int dimAc=(int)inputintec.n_cols;
    
    outputf  << "The length of the matrix is: " << dimA << endl;
    Amat.set_size(dimA,dimA);
    
    Amat.zeros();
    evawdlfn(c,0.0,wdlf1,atzero,negindex);
    
    int rows=(int)Amat.n_rows;
    
    int cols=(int)Amat.n_cols;
    
    int rowsv=(int)inputintec.n_cols;
    
    
    
    int chunk = int(floor(cols/OMP_NUM_THREADS));
    
    timer.tic();
    
#pragma omp parallel shared(Amat,inputintec,c,cols,wdlf1,wdlf2,rows,atzero,problem,normal,chunk) private(tid,j,k)
    {
        Row<double> diffsave(rowsv);
        Row<double> savingcallj(dimAc);
        Row<double> savingcallk(dimAc);
        Row<double> resultj(dimAc);
        Row<double> resultk(dimAc);
        
        diffsave.zeros();
        savingcallj.zeros();
        savingcallk.zeros();
        resultj.zeros();
        resultk.zeros();
        
        double twopointsdistance=0.0;
        
        double wdlfvalue1=0.0;
        double wdlfvalue2=0.0;
        double checking=0.0;

#if defined(_OPENMP)
        tid=omp_get_thread_num();
        if (tid == 0)
        {
            outputf << "Initializing matrices...\n" << endl;
        }
#endif
#pragma omp barrier
#pragma omp for schedule (dynamic, chunk) nowait
        for(j=0; j<cols; ++j)
        {
            savingcallj=inputintec(j, All);
            odesystem(problem,normal,savingcallj,resultj);
            for(k=0; k<rows; ++k)
            {
                savingcallk=inputintec(k, All);
                diffsave=savingcallj-savingcallk;
                if(k==j){
                    Amat(j,k)=-atzero*dot(resultj,resultj);
                }else{
                    twopointsdistance=sqrt(dot(diffsave,diffsave));
                    checking=1.0-c*twopointsdistance;
                    if(checking>0.0)
                    {
                        odesystem(problem,normal,savingcallk,resultk);
                        evawdlfn(c, twopointsdistance,wdlf1,wdlfvalue1,negindex);
                        
                        evawdlfn(c, twopointsdistance,wdlf2,wdlfvalue2,negindex);
                        
                        Amat(j,k)=-wdlfvalue2*dot(diffsave,resultj)*dot(diffsave,resultk)-wdlfvalue1*dot(resultj,resultk);
                    }
                }
            }
        }
    }
    outputf  << " " << endl;
    outputf  << "Computing the Interpolation matrix lasted: " << timer.toc() << endl;
    
    printhour(1);
}



void choldecom(mat &Amat, mat &R)
{
    wall_clock timer;
    timer.tic();
    
    

    int maxite=(int)Amat.n_rows;
    
    R.resize(maxite,maxite);
    R = chol( Amat );
    
    outputf  << " " << endl;
    outputf  << "The whole Cholesky proceedure last: " << timer.toc() << endl;
    
    printhour(1);
}

void orbitalderivative(glovar::Problema problem, int iteration, bool normal, double c, mat &R, Col<double> &alp, vec &betaod, Mat<double> &inputinterior, Row<double> &wedlf, Mat<double> &domaingrid, Row<double> &qmatod, Row<double> &qmatderod)
{
    
    

    bool negindex=false;
    
#if defined(_OPENMP)
    omp_set_dynamic(0);
    omp_set_num_threads(OMP_NUM_THREADS);
#endif
    int i=0, k=0, tid=0;
    
    wall_clock timer, timer1;
    
    timer.tic();
    if(normal)
    {
        outputf  << "Computing the Lyapunov function with almost normalized function" << endl;
    }else{
        outputf  << "Computing the Lyapunov function with no normalized function" << endl;
    }
    
    int chunk;
    int maxite=(int)domaingrid.n_rows;
    int maxbet=(int)inputinterior.n_rows;
    int pointdim=(int)inputinterior.n_cols;
    int wdfdim=(int)wedlf.size();
    
    
    qmatod.resize(maxite);
    qmatderod.resize(maxite);
    qmatod.zeros();//=0.0;
    qmatderod.zeros();//=0.0;
    
    Row<double> wdlf1(wdfdim-2), wdlf2(wdfdim-4);
    
    
    
    betaod.resize(maxbet);
    
    timer1.tic();
    betaod=solve(trimatu(R), solve(trimatl(trimatu(R).t()), alp),solve_opts::fast);
    outputf  << "The whole proceedure solve the Lyapunov equation lasted: " << timer1.toc() << endl;
    wendlandderivative(wedlf,wdlf1,negindex);
    wendlandderivative(wdlf1,wdlf2,negindex);
    
    chunk = int(floor(maxite/OMP_NUM_THREADS));
    
#pragma omp parallel num_threads(OMP_NUM_THREADS) shared(qmatod,qmatderod,chunk) private(tid,i,k)
    {//
#if defined(_OPENMP)
        tid=omp_get_thread_num();
        if(tid==0)
        {
            outputf  << "Starting the parallel region" << endl;
        }
#endif
        Row<double> diffpoints(pointdim), diffpointski(pointdim), diffpointskineg(pointdim);
        Row<double> resulti(pointdim), resultk(pointdim), saving(pointdim), savingdomain(pointdim);
        
        diffpoints.zeros();
        diffpointski.zeros();
        diffpointskineg.zeros();
        resulti.zeros();
        resultk.zeros();
        saving.zeros();
        savingdomain.zeros();
        
        double checking=0.0;
        double proctk=0.0;
        double producting=0.0;
        double twopointsdistance=0.0;
        double wdlfvalue1=0.0;
        double wdlfvalue2=0.0;
#if defined(_OPENMP)
        tid=omp_get_thread_num();
#endif
#pragma omp barrier
#pragma omp for schedule (dynamic, chunk) nowait
        for(i=0; i<maxite; ++i)
        {
            savingdomain=domaingrid(i,All);
            odesystem(problem,normal,savingdomain,resulti);
            for(k=0; k<maxbet; ++k)
            {
                saving=inputinterior(k,All);
                diffpoints=savingdomain-saving;
                twopointsdistance=sqrt(dot(diffpoints,diffpoints));
                checking=1.0-c*twopointsdistance;
                if(checking>0.0)
                {
                    odesystem(problem,normal,saving,resultk);
                    evawdlfn(c, twopointsdistance, wdlf1, wdlfvalue1,negindex);
                    evawdlfn(c, twopointsdistance, wdlf2, wdlfvalue2,negindex);
                    diffpointski=saving-savingdomain;
                    proctk=dot(diffpointski,resultk);
                    producting=betaod(k)*proctk;
                    qmatod(i)+=producting*wdlfvalue1;
                    qmatderod(i)+=-wdlfvalue2*producting*dot(diffpointski,resulti)-betaod(k)*wdlfvalue1*dot(resulti,resultk);
                }
            }
        }
    }
    if(glovar::printing)
    {
        printall("alpha", fextension, iteration, alp);
        
        printall("vecbetaod",fextension, iteration,betaod);
        
        if(normal)
        {
            printall("lyapdeltaite", fextension, iteration, qmatod);
            printall("lyapprimedeltaite", fextension, iteration, qmatderod);
        }else{
            printall("lyapite", fextension, iteration, qmatod);
            printall("lyapprimeite", fextension, iteration, qmatderod);
        }
        
    }
    
    outputf  << " " << endl;
    outputf  << "The whole proceedure to compute the Lyapunov function lasted: " << timer.toc() << endl;
    
    printhour(1);
}





