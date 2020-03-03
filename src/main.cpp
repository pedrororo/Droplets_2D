// Author: Lucas Berg

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "options/user_options.h"
#include "solver/droplets_solver.h"

using namespace std;

void configure_solver_from_options (struct solver_droplets *the_solver,\
									struct user_options *the_options)
{
	// Configure ODE solver
	double tmax = the_options->simulation_time;
	double dt_ode = the_options->dt_ode;
	configure_ode_solver(the_solver->the_ode_solver,tmax,dt_ode);
	
	// Configure PDE solver
	double dt_pde = the_options->dt_pde;
	double dx_pde = the_options->start_dx;
	double dy_pde = the_options->start_dy;
	double sigma_x = the_options->sigma_x;
	double sigma_y = the_options->sigma_y;
	uint32_t nx = nearbyint(the_options->side_length / the_options->start_dx);
	uint32_t ny = nearbyint(the_options->side_length / the_options->start_dy);
	uint32_t nt = nearbyint(the_options->simulation_time / the_options->dt_pde);
	configure_pde_solver(the_solver->the_pde_solver,sigma_x,sigma_y,dt_pde,dx_pde,dy_pde,nx,ny,nt);
	
	// Configure global solver
	the_solver->print_rate = the_options->print_rate;
	the_solver->spiral_time = the_options->spiral_time;
	if (the_options->output_dir != NULL)
	{
		uint32_t size = strlen(the_options->output_dir);
		the_solver->output_dir = (char*)malloc(sizeof(char)*size);
		strcpy(the_solver->output_dir,the_options->output_dir);
	}
	if (the_options->input_filename_U != NULL && the_options->input_filename_V != NULL)
	{
		uint32_t size;
		size = strlen(the_options->input_filename_U);
		the_solver->load_state_filename_U = (char*)malloc(sizeof(char)*size);
		strcpy(the_solver->load_state_filename_U,the_options->input_filename_U);
		
		size = strlen(the_options->input_filename_V);
		the_solver->load_state_filename_V = (char*)malloc(sizeof(char)*size);
		strcpy(the_solver->load_state_filename_V,the_options->input_filename_V);
		
		the_solver->using_load_state = true;
	}
	
}


int main (int argc, char *argv[])
{
    if (argc-1 != 1)
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    struct user_options *the_options = new_user_options(argc,argv);
    struct solver_droplets *the_solver = new_solver_droplets();
    
    configure_solver_from_options(the_solver,the_options);
    
    solve_problem(the_solver);

    //grow_tree(the_network,options);
    //write_to_vtk(the_network);

    //free_cco_network(the_network);
    free_user_options(the_options);

    return 0;
}
