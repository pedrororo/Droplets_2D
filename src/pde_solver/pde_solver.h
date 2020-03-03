
#ifndef PDE_SOLVER_H
#define PDE_SOLVER_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>

struct solver_pde
{
	double dtime;
	
	//LINEAR SYSTEM VARIABLES
	double *au_x; // inferior diagonal
	double *bu_x; // main diagonal
	double *cu_x; // superior diagonal
	double *av_x; // inferior diagonal
	double *bv_x; // main diagonal
	double *cv_x; // superior diagonal
	double *au_y; // inferior diagonal
	double *bu_y; // main diagonal
	double *cu_y; // superior diagonal
	double *av_y; // inferior diagonal
	double *bv_y; // main diagonal
	double *cv_y; // superior diagonal
	double *rhsu_x; // right hand side
	double *rhsu_y; // right hand side
	double *rhsv_x; // right hand side
	double *rhsv_y; // right hand side
	
	double *Aux_UX;
	double *Aux_UY;
	double *Aux_VX;
	double *Aux_VY;
	
	double  sigma_x;  
	double  sigma_y;  
	double  L0_x;  
	double  L0_y;  
	double  L1_x;  
	double  L1_y;

};

struct solver_pde* new_solver_pde();

#endif
