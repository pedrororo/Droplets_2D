#include "pde_solver.h"

struct solver_pde* new_solver_pde()
{
	
	struct solver_pde *pde_data = (struct solver_pde*)malloc(sizeof(struct solver_pde));
	
	return pde_data;
}

