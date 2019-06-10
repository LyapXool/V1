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


#include <list>
#include <armadillo>
#include "chainrecurrentsets.hpp"
#include "odesystem.hpp"
#include "RBF.hpp"
#include "instructions.hpp"
#include "wendland.hpp"
#include "generalities.hpp"

#define PI 3.141592653589793
using namespace arma;
using namespace std;
using namespace glovar;

arma::span const All=span::all;


void bigflowerh(const string gridtoeval,glovar::Problema problem, bool normal, int angles, int circles, bool spherical,  double maxmaxx, double minminx, double maxmaxy, double minminy, double maxmaxz, double minminz, double alpha, double radio, Mat<double> &RBFbasis, Mat<double> &inputinterior, Mat<double> &bigflowergrid, Mat<double> &biganalgrid, vector<bool> &boolbigflowergrid, vector<bool> &boolbiganalgrid, Mat<double> &cleanbigag, Mat<double> &cleanbigfg, int &stride)
{
    wall_clock timer;
    timer.tic();

    int lcols=(int)inputinterior.n_cols;
    int lrows=(int)inputinterior.n_rows;
    
    if(gridtoeval=="circular")
    {
        if(inputinterior.n_cols==2)
        {
            stride=circles*angles+1;
            
            bigflowergrid.resize(lrows*(stride-1),lcols);
            bigflowergrid.zeros();
            biganalgrid.resize(stride*lrows,lcols);
            biganalgrid.zeros();
            if(spherical)
            {
                for(int i=0; i<lrows; ++i)
                {
                    int j=stride*i;
                    biganalgrid(j,All)=inputinterior(i,All);
                    for(int q=0; q<circles; q++)
                    {
                        double div=((double)q+1.0)/(double)circles;
                        for(int k=0; k<angles;++k)
                        {
                            biganalgrid(j+k+1+q*angles,0)=div*radio*alpha*cos((2.0*PI*k)/angles)+inputinterior(i,0);
                            biganalgrid(j+k+1+q*angles,1)=div*radio*alpha*sin((2.0*PI*k)/angles)+inputinterior(i,1);
                        }
                    }
                }
                
                for(int i=0; i<lrows; ++i)
                {
                    
                    int j=(stride-1)*i;
                    for(int q=0; q<circles; q++)
                    {
                        double div=((double)q+1.0)/(double)circles;
                        for(int k=0; k<angles;++k)
                        {
                            bigflowergrid(j+k+q*angles,0)=div*radio*alpha*cos((2.0*PI*k)/angles)+inputinterior(i,0);
                            bigflowergrid(j+k+q*angles,1)=div*radio*alpha*sin((2.0*PI*k)/angles)+inputinterior(i,1);
                        }
                    }
                }

                
            }else{
                stride=2*lcols*angles+1;
                bigflowergrid.resize(2*lcols*angles*lrows,lcols);
                biganalgrid.resize(stride*lrows,lcols);
                double norm1=0.0;
                double norm2=0.0;
                bigflowergrid.zeros();
                biganalgrid.zeros();
                int j, jd;
                Mat<double> canonicalset(lcols,lcols);
                {
                    canonicalset=RBFbasis;
                }
                norm1=norm(canonicalset(All,0));
                norm2=norm(canonicalset(All,1));
                {
                    for(int i=0; i<lrows; ++i)
                    {
                        j=stride*i;
                        jd=(stride-1)*i;
                        biganalgrid(j,All)=inputinterior(i,All);
                        int kp=0;
                        for(int kd=0; kd<angles; kd+=1)
                        {
                            bigflowergrid(jd+kp,All)=inputinterior(i,All)+(radio/angles)*(kd+1)*alpha*(canonicalset(0,All)/norm1);
                            bigflowergrid(jd+kp+1,All)=inputinterior(i,All)-(radio/angles)*(kd+1)*alpha*(canonicalset(0,All)/norm1);
                            bigflowergrid(jd+kp+2,All)=inputinterior(i,All)+(radio/angles)*(kd+1)*alpha*(canonicalset(1,All)/norm2);
                            bigflowergrid(jd+kp+3,All)=inputinterior(i,All)-(radio/angles)*(kd+1)*alpha*(canonicalset(1,All)/norm2);

                            biganalgrid(j+kp+1,All)=bigflowergrid(jd+kp,All);
                            biganalgrid(j+kp+2,All)=bigflowergrid(jd+kp+1,All);
                            biganalgrid(j+kp+3,All)=bigflowergrid(jd+kp+2,All);
                            biganalgrid(j+kp+4,All)=bigflowergrid(jd+kp+3,All);

                            kp+=4;
                        }
                    }
                }
            }
        }
        if(inputinterior.n_cols==3){
            
            
            
            if(lcols!=3)
            {
                
                

                outputf  << "ERROR!! YOU ARE SUPPOSED TO BE RUNNING 3D-SYSTEMS" << endl;
                outputf  << "VERIFY THE DIMENSION OF YOUR PROBLEM!" << endl;
                
                
                exit(9);
            }
            Mat<double> transcoord(angles,lcols);
            bigflowergrid.zeros();
            biganalgrid.zeros();
            if(spherical)
            {
                bigflowergrid.resize(angles*lrows*circles,lcols);
                biganalgrid.resize((circles*angles+1)*lrows,lcols);
                stride=circles*angles+1;
                int j, jj;
                double phi, theta;
                {
                    for(int i=0; i<lrows; ++i)
                    {
                        j=stride*i;
                        jj=(stride-1)*i;
                        biganalgrid(j,All)=inputinterior(i,All);
                        for(int q=0; q<circles; q++)
                        {
                            double div=((double)q+1.0)/(double)circles;
                            for(int k=0; k<angles;++k)
                            {
                                phi=(acos(1.0-2.0*(k+0.5)/angles));
                                theta=(1+2.23606797749979*(k+0.5));
                                
                                transcoord(k,0)=div*radio*radio*alpha*cos(theta)*sin(phi);
                                transcoord(k,1)=div*radio*radio*alpha*sin(theta)*sin(phi);
                                transcoord(k,2)=div*radio*radio*alpha*cos(phi);
                                
                                bigflowergrid(jj+k+q*angles,All)=transcoord(k,All)+inputinterior(i,All);
                                biganalgrid(j+k+q*angles+1,All)=bigflowergrid(jj+k+q*angles,All);
                            }
                        }
                    }
                }
            }else{
                stride=2*lcols*angles+1;
                bigflowergrid.resize(2*lcols*angles*lrows,lcols);
                biganalgrid.resize(stride*lrows,lcols);
                double norm3=0.0;
                double norm1=0.0;
                double norm2=0.0;
                bigflowergrid.zeros();
                biganalgrid.zeros();
                int j, jd;
                Mat<double> canonicalset(lcols,lcols);
                {
                    canonicalset=RBFbasis/radio;
                }
                norm3=norm(canonicalset(All,0));
                norm1=norm(canonicalset(All,1));
                norm2=norm(canonicalset(All,2));
                {
                    for(int i=0; i<lrows; ++i)
                    {
                        j=stride*i;
                        jd=(stride-1)*i;
                        biganalgrid(j,All)=inputinterior(i,All);
                        int kp=0;
                        for(int kd=0; kd<angles; kd+=1)
                        {
                            bigflowergrid(jd+kp,All)=inputinterior(i,All)+(radio/angles)*(kd+1)*alpha*(canonicalset(0,All)/norm3);
                            bigflowergrid(jd+kp+1,All)=inputinterior(i,All)-(radio/angles)*(kd+1)*alpha*(canonicalset(0,All)/norm1);
                            biganalgrid(j+kp+1,All)=bigflowergrid(jd+kp,All);
                            biganalgrid(j+kp+2,All)=bigflowergrid(jd+kp+1,All);
                            bigflowergrid(jd+kp+2,All)=inputinterior(i,All)+(radio/angles)*(kd+1)*alpha*(canonicalset(1,All)/norm1);
                            bigflowergrid(jd+kp+3,All)=inputinterior(i,All)-(radio/angles)*(kd+1)*alpha*(canonicalset(1,All)/norm1);
                            biganalgrid(j+kp+3,All)=bigflowergrid(jd+kp+2,All);
                            biganalgrid(j+kp+4,All)=bigflowergrid(jd+kp+3,All);
                            bigflowergrid(jd+kp+4,All)=inputinterior(i,All)+(radio/angles)*(kd+1)*alpha*(canonicalset(2,All)/norm2);
                            bigflowergrid(jd+kp+5,All)=inputinterior(i,All)-(radio/angles)*(kd+1)*alpha*(canonicalset(2,All)/norm2);
                            biganalgrid(j+kp+5,All)=bigflowergrid(jd+kp+4,All);
                            biganalgrid(j+kp+6,All)=bigflowergrid(jd+kp+5,All);
                            kp+=6;
                        }
                    }
                }
            }
        }
    }
    if(gridtoeval=="directional")
    {
        int newlenght=(int)(angles*2*lrows);
        
        
        int j,jd;
        stride=angles*2+1;

        double norm;
        biganalgrid.resize(lrows*stride,lcols);
        bigflowergrid.resize(newlenght,lcols);
        
        
        Mat<double> domain(newlenght,lcols);
        Row<double> savingdomain(lcols), evaldfunction(lcols);
        {
            for(int i=0; i<lrows; ++i)
            {
                j=stride*i;
                jd=(stride-1)*i;
                biganalgrid(j,All)=inputinterior(i,All);
                savingdomain=inputinterior(i,All);
                
                odesystem(problem,normal,savingdomain, evaldfunction);
                norm=sqrt(dot(evaldfunction,evaldfunction));
                int kp=0;
                for(int kd=0; kd<angles; kd+=1)
                {
                    bigflowergrid(jd+kp,All)=inputinterior(i,All)+(radio/angles)*(kd+1)*alpha*(evaldfunction/norm);
                    bigflowergrid(jd+kp+1,All)=inputinterior(i,All)-(radio/angles)*(kd+1)*alpha*(evaldfunction/norm);
                    biganalgrid(j+kp+1,All)=bigflowergrid(jd+kp,All);
                    biganalgrid(j+kp+2,All)=bigflowergrid(jd+kp+1,All);
                    kp+=2;
                    
                }
            }
        }
     }
    
    list<int> counter,counterf;
    if(lcols==2)
    {
        for(int i=0; i<(int)biganalgrid.n_rows; ++i)
        {
            if((biganalgrid(i,0)<=maxmaxx && biganalgrid(i,0)>=minminx) && (biganalgrid(i,1)<=maxmaxy && biganalgrid(i,1)>=minminy))
            {
                counter.push_back(i);
            }
        }
        for(int ii=0; ii<(int)bigflowergrid.n_rows; ++ii)
        {
            if((bigflowergrid(ii,0)<=maxmaxx && bigflowergrid(ii,0)>=minminx) && (bigflowergrid(ii,1)<=maxmaxy && bigflowergrid(ii,1)>=minminy))
            {
                counterf.push_back(ii);
            }
        }
    }
    if(lcols==3)
    {
        for(int i=0; i<(int)biganalgrid.n_rows; ++i)
        {
            if((biganalgrid(i,0)<=maxmaxx && biganalgrid(i,0)>=minminx) && (biganalgrid(i,1)<=maxmaxy && biganalgrid(i,1)>=minminy) && (biganalgrid(i,2)<=maxmaxz && biganalgrid(i,2)>=minminz))
            {
                counter.push_back(i);
            }
        }
        for(int ii=0; ii<(int)bigflowergrid.n_rows; ++ii)
        {
            if((bigflowergrid(ii,0)<=maxmaxx && bigflowergrid(ii,0)>=minminx) && (bigflowergrid(ii,1)<=maxmaxy && bigflowergrid(ii,1)>=minminy) && (bigflowergrid(ii,2)<=maxmaxz && bigflowergrid(ii,2)>=minminz))
            {
                counterf.push_back(ii);
            }
        }
    }
    
    
    
    int ana=(int)counter.size();
    int flo=(int)counterf.size();
    
    boolbiganalgrid.resize(stride*lrows,false);
    boolbigflowergrid.resize(lrows*(stride-1),false);
    cleanbigag.resize(ana,lcols);
    cleanbigfg.resize(flo,lcols);
    
    int n=0;
    int m=0;
    
    for(list<int>::iterator i=counter.begin(); i!=counter.end(); ++i)
    {
        boolbiganalgrid[*i]=true;
        cleanbigag(n,All)=biganalgrid(*i,All);
        n++;
    }
    for(list<int>::iterator ii=counterf.begin(); ii!=counterf.end(); ++ii)
    {
        boolbigflowergrid[*ii]=true;
        cleanbigfg(m,All)=bigflowergrid(*ii,All);
        m++;
    }
    
    counter.clear();
    counterf.clear();
    
    
    
    
    if(glovar::printing)
    {
        int iteration=0;
        int fin=(int)(bigflowergrid.n_rows);
        
        if(lcols==2)
        {
            Col<double> xx(fin);
            Col<double> yy(fin);
            xx=bigflowergrid(All,0);
            yy=bigflowergrid(All,1);
            printall("x", fextension, iteration, xx);
            printall("y", fextension, iteration, yy);
        }else if(lcols==3)
        {
            Col<double> xx(fin);
            Col<double> yy(fin);
            Col<double> zz(fin);
            xx=bigflowergrid(All,0);
            yy=bigflowergrid(All,1);
            zz=bigflowergrid(All,2);
            printall("x", fextension, iteration, xx);
            printall("y", fextension, iteration, yy);
            printall("z", fextension, iteration, zz);
        }
    }
    
    
    
    
    outputf  << "The whole proceedure to construct the evaluation grid lasted: " << timer.toc() << endl;
    
    
}






void qpzerpctotalh(int iteration, int angles, int circles, int &stride , const int &defcase, double critval, Row<double> &qfunc, Row<double> &qmatderod, Row<double> &failingq, Row<double> &failingqprime,  Mat<double> &angrid, Mat<double> &inputinterior, Col<double> &newalpha, Mat<double> &failinggrid, Mat<double> &angridnonzero, Mat<double> &flowernozero, Mat<double> &collocationzero, Mat<double> &collocationnozero, Row<double> &alphazero, Row<double> &alphanozero,list<int> &counter, list<int> &counterncol,vector<bool> &boolbigflowergrid)
{
    wall_clock timer;
    timer.tic();
    
    
    
    counterncol.clear();
    counter.clear();
    
    if(stride<=0.5)
    {
        outputf  << " Relation between the evaluation grid and the collocation grid has failed. Obtaining the correct relation in a moment and continuing... " << endl;
        stride=circles*angles+1;
    }
    list<int> countern;
    countern.clear();
    
    int maxite=(int)qmatderod.n_cols;
    int pointdim=(int)angrid.n_cols;
    
    
    int div, divc,divco;
    divco=-1;
    int  i, k;
    {
        for(i=0; i<maxite; i+=(stride-1))
        {
            
            bool checking=true;
            for(k=0; k<(stride-1); ++k)
            {
                if(qmatderod(i+k)>critval && boolbigflowergrid[i+k]){
                    checking=false;
                    countern.push_back(i+k);
                    divc=trunc(i/(stride-1));
                    if(divco!=divc){
                        counterncol.push_back(divc);
                    }
                    divco=divc;
                }
            }
            
            if( checking )
            {
                div=trunc(i/(stride-1));
                counter.push_back(div);
            }
        }
    }
    
    int counternsize=(int)countern.size();
    int countersize=(int)counter.size();
    int alphasize=(int)trunc((int)angrid.n_rows/stride);
    int colzerosize=(int)counterncol.size();//int colzerosize=(int)abs(alphasize-countersize);
    failinggrid.resize(counternsize,pointdim);
    failingq.resize(counternsize);
    failingqprime.resize(counternsize);
    angridnonzero.resize(stride*countersize,pointdim);
    flowernozero.resize((stride-1)*countersize,pointdim);
    collocationnozero.resize(countersize,pointdim);
    alphanozero.resize(countersize);
    collocationzero.resize(colzerosize,pointdim);
    alphazero.resize(colzerosize);
    
    
    
    if(counternsize==0)
    {
        outputf  << "Under the given settings and critical value, the algorithm has not found any point in the chain-recurrent set. Please, restart the calculation with a stricter critical value" << endl;
        exit(9);
    }
    int n=0;
    int m=0;
    int p=0;
    
    for(list<int>::iterator i=counter.begin(); i!=counter.end(); ++i)
    {
        for(int j=0; j<=(stride-1);++j)
        {
            angridnonzero(stride*n+j,All)=angrid(stride*(*i)+j,All);
        }
        
        for(int k=1; k<=(stride-1);++k)
        {
            flowernozero((stride-1)*n+k-1,All)=angrid(stride*(*i)+k,All); //flowernozero comienza en cero cuando angrid comienza en 1
        }
        
        collocationnozero(n,All)=angrid(stride*(*i),All);
        
        n++;
    }
   
    for(list<int>::iterator ii=countern.begin(); ii!=countern.end(); ++ii)
    {
        failinggrid(m,All)=inputinterior((*ii), All);
        failingq(m)=qfunc((*ii));
        failingqprime(m)=qmatderod((*ii));
        m++;
    }
    
    for(list<int>::iterator iii=counterncol.begin(); iii!=counterncol.end(); ++iii)
    {
        collocationzero(p,All)=angrid(stride*(*iii),All);
        p++;
    }
    
    if(defcase==1)
    {
        newalpha.resize(alphasize);
        newalpha.fill(-1.0);
        int mm=0;
        for(list<int>::iterator ii=counterncol.begin(); ii!=counterncol.end(); ++ii)
        {
            newalpha((*ii))=0.0;
            alphazero(mm)=0.0;
            mm++;
        }
        alphanozero.fill(-1.0);
    }
    
    if(defcase==2)
    {
        int vlengthnz=(int)collocationnozero.n_rows;
        int vlengthz=(int)collocationzero.n_rows;
        Row<double> saving(pointdim);
        Row<double> distances(vlengthnz);
        distances.fill(1000000.0000);
        double d;
        for(int jj=0; jj<vlengthnz; ++jj)
        {
            for(int j=0; j<vlengthz; ++j)
            {
                saving=collocationzero(j,All)-collocationnozero(jj,All);
                d=sqrt(dot(saving,saving));
                if(d<distances(jj))
                {
                    distances(jj)=d;
                }
            }
        }
        
        newalpha.resize(alphasize);
        for(list<int>::iterator iii=counterncol.begin(); iii!=counterncol.end(); ++iii)
        {
            newalpha(*iii)=0.0;
        }
        alphazero.zeros();
        int n=0;
        for(list<int>::iterator i=counter.begin(); i!=counter.end(); ++i)
        {
            newalpha(*i)=-exp(-1.0/(300.0*distances(n)*distances(n)));
            alphanozero(n)=newalpha(*i);
            n++;
        }
    }
    
    if(defcase==3)
    {
        double sumgood;
        double sumbads;
        newalpha.resize(alphasize);
        int pa=0;
        for(list<int>::iterator iii=counterncol.begin(); iii!=counterncol.end(); ++iii)
        {
            sumbads=0.0;
            for(int j=0; j<(stride-1);++j)
            {
                sumbads+=qmatderod((stride-1)*(*iii)+j);
                
            }
            if(sumbads>0.0)
            {
                sumbads=0.0;
            }
            newalpha(*iii)=sumbads/((double)(stride-1));
            alphazero(pa)=newalpha(*iii);
            pa++;
        }
        int na=0;
        for(list<int>::iterator i=counter.begin(); i!=counter.end(); ++i)
        {
            sumgood=0.0;
            for(int j=0; j<(stride-1);++j)
            {
                sumgood+=qmatderod((stride-1)*(*i)+j);
            }
            if(sumgood>0.0)
            {
                sumgood=0.0;
            }
            newalpha(*i)=sumgood/((double)(stride-1));
            alphanozero(na)=newalpha(*i);
            na++;
        }
        
        double normalizationfactor=0.0;
        for(int i=0; i<alphasize; ++i)
        {
            normalizationfactor+=newalpha(i);
        }
        
        newalpha=abs(alphasize/normalizationfactor)*newalpha;
        double prueba=0.0;
        for(int i=0; i<alphasize; ++i)
        {
            prueba+=newalpha(i);
        }
    }
    
    
    counterncol.clear();
    counter.clear();
    countern.clear();
    
    
    
    if(glovar::printing)
    {
        if(pointdim==2)
        {
            
            Col<double> x(counternsize);
            Col<double> y(counternsize);
            Col<double> z(counternsize);
            Row<double> w(counternsize);
            Row<double> v(counternsize);
            Col<double> t(countersize);
            Col<double> s(countersize);
            Col<double> r(countersize);
            
            
            Col<double> u(colzerosize);
            Col<double> q(colzerosize);
            Col<double> o(colzerosize);
            Col<double> l(colzerosize);
            Col<double> h(colzerosize);
            
            x=failinggrid(All,0);
            y=failinggrid(All,1);
            
            w=failingq;
            
            v=failingqprime;
            t=collocationnozero(All,0);
            s=collocationnozero(All,1);
            
            u=collocationzero(All,0);
            q=collocationzero(All,1);
            
            printall("fx", fextension, iteration, x);
            printall("fy", fextension, iteration, y);
            printall("flf", fextension, iteration, w);
            printall("flfp", fextension, iteration, v);
            printall("xcollocationnozero", fextension, iteration, t);
            printall("ycollocationnozero", fextension, iteration, s);
            printall("xcollocationzero", fextension, iteration, u);
            printall("ycollocationzero", fextension, iteration, q);
        }
       
        if(pointdim==3)
        {
            Col<double> x(counternsize);
            Col<double> y(counternsize);
            Col<double> z(counternsize);
            Row<double> w(counternsize);
            Row<double> v(counternsize);
            Col<double> t(colzerosize);
            Col<double> s(colzerosize);
            Col<double> r(colzerosize);
            Col<double> u(colzerosize);
            Col<double> q(colzerosize);
            Col<double> o(colzerosize);
            Col<double> l(colzerosize);
            Col<double> h(colzerosize);
            
            
            
            x=failinggrid(All,0);
            y=failinggrid(All,1);
            z=failinggrid(All,2);
            w=failingq;
            v=failingqprime;

            
            t=collocationnozero(All,0);
            s=collocationnozero(All,1);
            r=collocationnozero(All,2);
            u=collocationzero(All,0);
            q=collocationzero(All,1);
            o=collocationzero(All,2);
            
            printall("fx", fextension, iteration, x);
            printall("fy", fextension, iteration, y);
            printall("fz", fextension, iteration, z);
            printall("flf", fextension, iteration, w);
            printall("flfp", fextension, iteration, v);
            
            printall("xcollocationnozero", fextension, iteration, t);
            printall("ycollocationnozero", fextension, iteration, s);
            printall("zcollocationnozero", fextension, iteration, r);
            printall("xcollocationzero", fextension, iteration, u);
            printall("ycollocationzero", fextension, iteration, q);
            printall("zcollocationzero", fextension, iteration, o);
        }
    }
    outputf  << " " << endl;
    outputf  << "The whole proceedure to filter the values lasted: " << timer.toc() << endl;
    
    printhour(1);
    
}
