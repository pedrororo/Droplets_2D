#include "ode_solver.h"

struct solver_ode* new_solver_ode()
{
	
	struct solver_ode *ode_data = (struct solver_ode*)malloc(sizeof(struct solver_ode));
	
	return ode_data;
	
	
}

void configure_ode_solver (struct solver_ode *s, const double tmax, const double dt)
{
	s->tmax = tmax;
	s->dt = dt;
}

void solve_ode_system_FHN (double **U, double **V, const double dt, const uint32_t nx, const uint32_t ny)
{
	double Aux_u, Aux_v, Iion_u, Iion_v;
	
	for (uint32_t i = 0; i < ny; i++)
	{
		for (uint32_t j = 0; j < nx; j++)
		{
			Iion_u = k_F*(U[i][j]*(1.0 - U[i][j])*(U[i][j] - a) - U[i][j]*V[i][j]); 
			Iion_v = k_F*test_K*(b*U[i][j] - V[i][j]);
			
			Aux_u  = (dt)*Iion_u;	
			Aux_v  = (dt)*Iion_v;
			
			U[i][j] = Aux_u + U[i][j];
			V[i][j] = Aux_v + V[i][j];
		}
	}
}
