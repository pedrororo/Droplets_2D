
#ifndef PDE_SOLVER_H
#define PDE_SOLVER_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>

struct solver_pde
{
	uint32_t nt;   //total iterations time
	
	double dt_pde;
	
	double dx_pde;
	double dy_pde;
	
	
	double sigma_x; //conductivity x
	double sigma_y; //conductivity y
	
	uint32_t    nx;		// Number of cells in grid x
	uint32_t    ny;		// Number of cells in grid y
	
	double L_x;  	//length x
	double L_y;     //length y
	
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
	 
	double  L0_x;  
	double  L0_y;  
	double  L1_x;  
	double  L1_y;

};

struct solver_pde* new_solver_pde();

void configure_pde_solver(struct solver_pde *s, const double sigma_x, const double sigma_y,\
						const double dt_pde, const double dx_pde, const double dy_pde,const uint32_t nx, const uint32_t ny, const uint32_t nt);

void allocate_pde_vectors (struct solver_pde *s);

void calculate_parabolic_LHS (struct solver_pde *s);

void set_parabolic_matrix (double *a, double *b, double *c, const double sigma, const double dt, const double h,const uint32_t n);

void solve_pde (struct solver_pde *the_pde_solver, double **U, double **V);

void solve_parabolic (const double *a, const double *b, const double *c, double *rhs, double **A, const uint32_t nx, const uint32_t ny, const char axis);

void solve_tridiagonal_system (const double *a, const double *b, const double *c, double *rhs, double *x, const uint32_t n);

#endif
