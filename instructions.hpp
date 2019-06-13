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
 ->
 */


#ifndef instructions_hpp
#define instructions_hpp


#include <stdio.h>
#include <iostream>

namespace glovar {
    extern int functionodecalls;

    enum Problema {
        TWOORBITS,VDP,HOMOCLINIC,DECREASING,TD1,TD2,TD3,SIMPLE3D,LORENZ
    };

    char const probnames[][11]={"TWOORBITS","VDP","HOMOCLINIC","DECREASING","TD1","TD2","TD3","SIMPLE3D","LORENZ"};


    /*%%%% SECTION TO DEFINE THE RBF AND EVALUATION GRIDS %%%%*/

    const int maxnegative=-450; /*VALUES FOR THE REGULAR CARTESIAN GRID TO BE TRANSFORMED*/
    const int maxpositive=450;

    const bool cartesian=false; /*WILL THE BASIS BE CANONICA=true OR HEXAGONAL=false*/

    const double alpha=0.1;//0.0665;//0165;  /*%% FOR HEXAGONAL BASIS, VALUE OF ALPHA. FOR CANONICAL BASIS WILL NOT BE USED.*/

    const double maxmaxx=1.0;   /* DEFINITION OF THE RBF COLLOCATION GRID BOUNDARIES*/
    const double minminx=-maxmaxx;
    const double maxmaxy=1.0;
    const double minminy=-maxmaxy;
    const double maxmaxz=1.0;
    const double minminz=-maxmaxz;

    const std::string gridtoeval="directional";  /*DIRECTIONAL OR CIRCULAR/SPHERICAL GRID*/

    const int circles=2; /*CIRFUMFERENCE FOR CIRCULAR ONLY*/
    const int angles=10; /* AMOUNT OF ANGLES FOR CIRCULAR OR POINTS FOR DIRECTIONAL*/

    const bool spherical=false; /*If it does not use the directional grid, then will it be spherical? False means along the normal axis set*/

    /*%%%% SECTION TO DEFINE PROBLEM TO ANALYZE %%%%*/

    const Problema problem=glovar::LORENZ;

    const int variable=3;

    const int ncritical=3;


    /*%%%% SECTION TO DEFINE CONDITIONS %%%%*/

    const bool constante=true; /*true IF ALPHA WILL BE -1 AT EVERY POINT AT THE BEGINNING. false if it will be -\|f(x)\|^2?*/

    const bool normal=true;    /*true FOR THE ALMOST NORMALIZED METHOD*/

    const int defcase=3;       /* REITERATION OF ALPHA:  1.- 0 and -1, 2.-exponential, 3.-average;*/


    const bool eigenvaluesjudge=false;                        /* analyzes the critical points*/

    const bool printing=true;


    const double critval=-0.5;

    const double radio=0.49;


    /*%%%% SECTION TO DEFINE AMOUNT OF ITERATIONS %%%%*/

    const int totaliterations=3;                             /*If it will iterate, how many will give?*/


    /*%%%% SECTION TO DEFINE WENDLAND FUNCTION %%%%*/

    const double l=4;
    const double k=2;
    const double c=1;


    /*%%%% GENERALITIES %%%%*/

    const int OMP_NUM_THREADS=16;

    const std::string fextension="m";


    extern std::ofstream outputf;

};




#endif /* instructions_hpp */
