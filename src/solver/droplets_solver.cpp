#include "droplets_solver.h"

struct solver_droplets* new_solver_droplets()
{
	struct solver_droplets *result = (struct solver_droplets*)malloc(sizeof(struct solver_droplets));
	
	result->the_ode_solver = new_solver_ode();
	result->the_pde_solver = new_solver_pde();
	
	result->output_dir = NULL;
	result->save_state_filename = NULL;
	result->load_state_filename_U = NULL;
	result->load_state_filename_V = NULL;
	result->using_load_state = false;
	
	return result;
}

void solve_problem (struct solver_droplets *s)
{
	allocate_common_vectors(s);
	allocate_pde_vectors(s->the_pde_solver);
	
	if (s->using_load_state == true)
	{
		printf("Reading from file '%s' and '%s'\n",s->load_state_filename_U,s->load_state_filename_V);
		read_steady_state_from_file(s);
	}
	else
	{
		printf("Default initial conditions\n");
		set_initial_conditions_default(s);
	}
	
	// TODO: Refine grid will enter here ...
	
	// Set the filename for the trace
	char output_trace_filename[200];
	sprintf(output_trace_filename,"%s/UU.dat",s->output_dir);

	// Fill the discretization matrix
	calculate_parabolic_LHS(s->the_pde_solver);
	
	double dt_pde = s->the_pde_solver->dt_pde;
	double dt_ode = s->the_ode_solver->dt;
	uint32_t print_rate = s->print_rate;
	uint32_t num_iterations_solve_ode = nearbyint(dt_pde/dt_ode);
	
	// TODO:
	// Set trace cells position
	//	**************
	//	* 	*	 *   *
	//    1		   2
	//	*	*    *   *
	//	**************
	//	*	*	 *	 *
	//	*	*	 *	 *
	//	**************
	//	* 	*	 *	 *
	//	  0		   3
	//	* 	*	 *   *
	//	**************
	
	// MAIN TIME LOOP
	bool first_call = true;
	uint32_t nx = s->the_pde_solver->nx;
	uint32_t ny = s->the_pde_solver->ny;
	uint32_t nt = s->the_pde_solver->nt;
	double dx = s->the_pde_solver->dx_pde;
	double dy = s->the_pde_solver->dy_pde;
	double **UU = s->UU;
	double **VV = s->VV;
	double spiral_time = s->spiral_time;
	char filename[200];
	char *output_dir = s->output_dir;
	for (uint32_t k = 0; k < nt; k++)
	{
		double t = k*dt_pde;
		
		for (uint32_t i = 0; i < num_iterations_solve_ode; i++)
		{
			solve_ode_system_FHN(UU,VV,dt_ode,nx,ny);
		}
		
		solve_pde(s->the_pde_solver,UU,VV);
		
		
		if (t >= spiral_time && first_call && s->using_load_state == false)
		{
			set_spiral(UU,dx,dy,nx,ny);
			first_call = false;
		}
		
		
		if (k % print_rate == 0)
		{
			printf("Time = %g\n",t);
			
			// Write the activator data
			sprintf(filename,"%s/UU_%d.dat",output_dir,k);
			write_result_to_file(filename,UU,nx,ny);
			
			// Write the inhibitor data
			sprintf(filename,"%s/VV_%d.dat",output_dir,k);
			write_result_to_file(filename,VV,nx,ny);
		}
	}
	

}

void allocate_common_vectors (struct solver_droplets *s)
{
	uint32_t nx = s->the_pde_solver->nx;
	uint32_t ny = s->the_pde_solver->ny;
	
	s->UU = (double**)malloc(sizeof(double*)*ny);
	s->VV = (double**)malloc(sizeof(double*)*ny);
	//s->UU_b = (double**)malloc(sizeof(double*)*ny);
	//s->VV_b = (double**)malloc(sizeof(double*)*ny);
	for (uint32_t i = 0; i < ny; i++)
	{
		s->UU[i] = (double*)malloc(sizeof(double)*nx);
		s->VV[i] = (double*)malloc(sizeof(double)*nx);
	}
}

void set_initial_conditions_default(struct solver_droplets *s)
{
	uint32_t nx = s->the_pde_solver->nx;
	uint32_t ny = s->the_pde_solver->ny;
	uint32_t nx_limit = nx / 10;
	double **UU = s->UU;
	double **VV = s->VV;
	
	for (uint32_t i = 0; i < ny; i++)
	{
		for (uint32_t j = 0; j < nx; j++)
		{
			if (j < nx_limit)
			{
				UU[i][j] = 0.9;
				VV[i][j] = 0.0;
			}
			//else
			//{
			//	UU[i][j] = 0.0;
			//	VV[i][j] = 0.0;
			//}
		}
		
	}
}

void set_spiral (double **U, const double dx, const double dy, const uint32_t nx, const uint32_t ny)
{
	for (uint32_t i = 0; i < ny; i++)
	{
		double y = i*dy;
		for (uint32_t j = 0; j < nx; j++)
		{
			double x = j*dx;
			if ((x >= 0.0 && x <= 0.0100) && (y >= 0.0 && y <= 0.0100))
			//~ if ((x >= 0.0 && x <= 0.0125) && (y >= 0.0 && y <= 0.0125))
			//~ if ((x >= 0.0 && x <= 0.0140) && (y >= 0.0 && y <= 0.0140))
			{
				U[i][j] = 0.9;
			}
		}
		
	}
}

void read_steady_state_from_file (struct solver_droplets *s)
{
	uint32_t nx = s->the_pde_solver->nx;
	uint32_t ny = s->the_pde_solver->ny;
	double **UU = s->UU;
	double **VV = s->VV;
	
	char *load_state_filename_U = s->load_state_filename_U;
	char *load_state_filename_V = s->load_state_filename_V;
	
	double value;
	
	FILE *file_U = fopen(load_state_filename_U,"r");
	FILE *file_V = fopen(load_state_filename_V,"r");
	for (uint32_t i = 0; i < ny; i++)
	{
		for (uint32_t j = 0; j < nx; j++)
		{
			fscanf(file_U,"%lf",&value);
			UU[i][j] = value;
			
			fscanf(file_V,"%lf",&value);
			VV[i][j] = value;
		}
	}
	fclose(file_U);
	fclose(file_V);
}

void write_result_to_file (const char filename[], double **A, const uint32_t nx, const uint32_t ny)
{
	FILE *file = fopen(filename,"w+");	
	
	for (uint32_t i = 0; i < ny; i++)	
	{
		for (uint32_t j = 0; j < nx; j++)
		{	
			fprintf(file, "%e\t", A[i][j]);	
		}
		fprintf(file,"\n"); fflush(file);
	}
	
	fclose(file);
}
