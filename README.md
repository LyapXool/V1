# V1
First version of LyapXool
To make this program work just do the next:
1) edit the odesystem.cpp and add the dynamical system you want to, using c++ language syntax. 
2) name it accordingly and add that name in instructions.hpp in the enum Problema { and in the char const probnames[][11]={
Follow the given sintax.
3) prepare the Makefile with the right libraries
4) edit the instructions.hpp file and modify:
     
     To use a hexagonal grid, keep false, to use cartesian, use true:
    const bool cartesian=false; 
     
     To vary the amount the collocation points, change the alpha parameter. The higher the value the lower the amount of collocation points:
    const double alpha=0.05;
     
     To define the phase space, remember that you can use independent values for the maximum and minimum in x,y and z:
    const double maxmaxx=1.5;   
    const double minminx=-maxmaxx;
    const double maxmaxy=1.5;
    const double minminy=-maxmaxy;
    const double maxmaxz=1.5;
    const double minminz=-maxmaxz;

     To define the type of evaluation grid. It can be directional, circular or spherical if it is a 3D case:
    const std::string gridtoeval="circular";  

    To modify the amount of concentric circles, please change the next value.
    const int circles=3; /*CIRFUMFERENCE FOR CIRCULAR ONLY*/
    
    To modify the amount of angles please change the next value. For the directional remember then that it will define the amount of points per direction:
    const int angles=100; 

    If it does not use the directional grid, then will it be spherical? False means along the normal axis set:
    const bool spherical=false;

    

    Select the problem:
    const Problema problem=glovar::TWOORBITS;

    Set the right dimension:
    const int variable=2;

    Set the amounf of critical points:
    const int ncritical=3;


   
    Chose if the firt alpha guess will be all equal to -1:
    const bool constante=false; 

    If you want to normalize the right handside of the ODE:
    const bool normal=false;    

    For reiterations: with 1 and 0, use 1, for exponential decay use 2 and for averaging use 3:
    const int defcase=3;      



    Do you want to analyse the critical points? Remember to defined the critical points in problem.cpp, line 75:
    const bool eigenvaluesjudge=false;                        /* analyzes the critical points*/

    Do you want the results printed? They will get printed in a file:
    const bool printing=true;

    What is the tolerance parameter for the chain-recurrent set:
    const double critval=-0.5;

    What is the radio for the evaluation grid?
    const double radio=0.499;



    How many iterations do you want?
    const int totaliterations=20;                             /*If it will iterate, how many will give?*/


   Conditions for the Wendland functions:

    const double l=5;
    const double k=3;
    const double c=1;


    
    Amount of processors:
    const int OMP_NUM_THREADS=16;

    Extension of the pritting file:
    const std::string fextension="m";


    extern std::ofstream outputf;

5) make and run
