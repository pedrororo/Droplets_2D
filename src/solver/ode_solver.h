
#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>

// FHN CONSTANT
const double a = 0.2;
const double b = 0.5;
const double k_F = 80.0;
const double test_K = 0.0011;

struct solver_ode
{
	
	double dt; 	
	double tmax;
	double time_new; 

};

struct solver_ode* new_solver_ode();

void configure_ode_solver (struct solver_ode *s, const double tmax, const double dt);

void solve_ode_system_FHN (double **U, double **V, const double dt, const uint32_t nx, const uint32_t ny);

#endif
