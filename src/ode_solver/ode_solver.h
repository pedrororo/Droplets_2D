
#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>

struct solver_ode
{
	
	double time; 	
	double dtime;
	double time_new;
	double a;
	double b;
	double **UU;
	double **VV;
	double **UU_b;
	double **VV_b;
	double sigma_x; //conductivity x
	double sigma_y; //conductivity y
	double L_x;  	//length x
	double L_y;     //length y
	int    sizet;   //total iterations time
	int    size_x;
	int    size_y; 

};

struct solver_ode* new_solver_ode();

void solveODESystem_FHN(struct solver_ode *s);

#endif
