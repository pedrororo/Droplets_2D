
#ifndef DROPLETS_SOLVER_H
#define DROPLETS_SOLVER_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <cmath>


#include "ode_solver.h"
#include "pde_solver.h"

struct solver_droplets
{
	double **UU;			// Matrix for the activator
	double **VV;			// Matrix for the inhibitor
	double **UU_b;			// Matrix for the activation (refined)
	double **VV_b;			// Matrix for the inhibitor (refined)
	
	uint32_t save_rate;
	char *save_state_filename;
	
	bool using_load_state;
	char *load_state_filename_U;
	char *load_state_filename_V;

	bool using_interpolation;	
	char *method_name;
	
	char *output_dir;
	uint32_t print_rate;
	
	double spiral_time;
	
	struct solver_ode *the_ode_solver;
	struct solver_pde *the_pde_solver;

};

struct solver_droplets* new_solver_droplets();

void solve_problem (struct solver_droplets *s);

void set_initial_conditions_default(struct solver_droplets *s);

double mh4(const double s1,const double s2,const double s3,const double s4);

void count_rows_and_columns (const char filename[], uint32_t &num_cell_refined_nx, uint32_t &num_cell_refined_ny);

void allocate_common_vectors (struct solver_droplets *s);

void set_spiral (double **U, const double dx, const double dy, const uint32_t nx, const uint32_t ny);

void read_steady_state_from_file (struct solver_droplets *s);

void write_result_to_file (const char filename[], double **A, const uint32_t nx, const uint32_t ny);

void write_result_interpolation_to_file (const char filename[], double **A, const double dx, const double dy, const uint32_t nx, const uint32_t ny);

void write_result_to_vtk (const char filename[], double **AA, const double dx, const double dy, const uint32_t nx, const uint32_t ny);

#endif
