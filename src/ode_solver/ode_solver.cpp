#include "ode_solver.h"

struct solver_ode* new_solver_ode()
{
	
	struct solver_ode *ode_data = (struct solver_ode*)malloc(sizeof(struct solver_ode));
	
	return ode_data;
	
	
}

