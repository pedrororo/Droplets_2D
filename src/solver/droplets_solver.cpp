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
	result->method_name = NULL;
	result->using_interpolation = false;
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
	char filename[200], filename_vtk[200];
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
			sprintf(filename_vtk,"%s/vtk/UU_%d.vtk",output_dir,k);
			if (s->using_interpolation == false)
				write_result_to_file(filename,UU,nx,ny);
			else
				write_result_interpolation_to_file(filename,UU,dx,dy,nx,ny);
			write_result_to_vtk(filename_vtk,UU,dx,dy,nx,ny);
			
			// Write the inhibitor data
			sprintf(filename,"%s/VV_%d.dat",output_dir,k);
			sprintf(filename_vtk,"%s/vtk/VV_%d.vtk",output_dir,k);
			if (s->using_interpolation == false)
				write_result_to_file(filename,VV,nx,ny);
			else
				write_result_interpolation_to_file(filename,VV,dx,dy,nx,ny);
			write_result_to_vtk(filename_vtk,VV,dx,dy,nx,ny);
		}
	}
	

}

void allocate_common_vectors (struct solver_droplets *s)
{
	uint32_t nx = s->the_pde_solver->nx;
	uint32_t ny = s->the_pde_solver->ny;
	
	s->UU = (double**)malloc(sizeof(double*)*ny);
	s->VV = (double**)malloc(sizeof(double*)*ny);
	s->UU_b = NULL;
	s->VV_b = NULL;
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
	double **UU_refined = s->UU_b;
	double **VV_refined = s->VV_b;
	double aux_U, aux_V;
	
	
	char *load_state_filename_U = s->load_state_filename_U;
	char *load_state_filename_V = s->load_state_filename_V;
	
	bool using_interpolation = s->using_interpolation;
	
	double value;
	char *aux_value;
	
	FILE *file_U;
	FILE *file_V;
	
	
	if(using_interpolation == true)
	{
		printf("[!] Using interpolation !!!\n");
		
		// Read the number of cells in each direction from the refined mesh
		uint32_t num_cell_refined_nx, num_cell_refined_ny;
		count_rows_and_columns(load_state_filename_U,num_cell_refined_nx,num_cell_refined_ny);
		
		// Allocate memory for the refined mesh
		UU_refined = (double**)malloc(sizeof(double*)*num_cell_refined_ny);
		VV_refined = (double**)malloc(sizeof(double*)*num_cell_refined_ny);
		for (uint32_t i = 0; i < num_cell_refined_ny; i++)
		{
			UU_refined[i] = (double*)malloc(sizeof(double)*num_cell_refined_nx);
			VV_refined[i] = (double*)malloc(sizeof(double)*num_cell_refined_nx);
		}
		
		// Read the values from the refined mesh
		file_U = fopen(load_state_filename_U,"r");
		file_V = fopen(load_state_filename_V,"r");
		for (uint32_t i = 0; i < num_cell_refined_ny; i++)
		{
			for (uint32_t j = 0; j < num_cell_refined_nx; j++)
			{
				fscanf(file_U,"%lf",&value);
				UU_refined[i][j] = value;
				
				fscanf(file_V,"%lf",&value);
				VV_refined[i][j] = value;
			}
		}
		fclose(file_U);
		fclose(file_V);
					
		//~ uint32_t ratio_x = (num_cell_refined_nx / nx);
		//~ uint32_t ratio_y = (num_cell_refined_ny / ny);
		double ratio_x = ((double)num_cell_refined_nx / (double)nx);
		double ratio_y = ((double)num_cell_refined_ny / (double)ny);
		
		// DEBUG
		printf("[interpolation] num_cell_refined_nx = %d\n",num_cell_refined_nx);
		printf("[interpolation] num_cell_refined_ny = %d\n",num_cell_refined_ny);
		printf("[interpolation] nx = %d\n",nx);
		printf("[interpolation] ny = %d\n",ny);
		printf("[interpolation] ratio_x = %f\n",ratio_x);
		printf("[interpolation] ratio_y = %f\n",ratio_y);
					
		for(uint32_t i = 0; i < ny-1; i++) 
		{
			uint32_t ib = i*ratio_y;
			for(uint32_t j = 0; j < nx-1; j++) 
			{
				uint32_t jb = j*ratio_x;		
				
				UU[i][j] = (mh4(UU_refined[ib][jb],UU_refined[ib][jb+1],UU_refined[ib+1][jb],UU_refined[ib+1][jb+1]));
				VV[i][j] = (mh4(VV_refined[ib][jb],VV_refined[ib][jb+1],VV_refined[ib+1][jb],VV_refined[ib+1][jb+1]));	
			}

		}	
	}	
	else
	{
		
		file_U = fopen(load_state_filename_U,"r");
		file_V = fopen(load_state_filename_V,"r");
		
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

}

double mh4(const double s1,const double s2,const double s3,const double s4)
{
	double result;
	
	//~ result  = (4.0*(s1*s2*s3*s4))/(s1+s2+s3+s4);
	result  = (4.0*(s1*s2*s3*s4))/(s1*s2*s3+s2*s3*s4+s3*s4*s1+s4*s1*s2);
	//~ result  = (2.0*(s1*s2))/(s1+s2);
	
	return result;
}

void count_rows_and_columns (const char filename[], uint32_t &num_cell_refined_nx, uint32_t &num_cell_refined_ny)
{
	FILE *file_U = fopen(filename,"r");
	
	int i;
	num_cell_refined_nx = 0;
	num_cell_refined_ny = 0;
	bool first_call = true;
	
	while((i = fgetc(file_U)) != EOF)
	{
		if((i=='\t')&&(first_call != false))
		{
			num_cell_refined_nx++;
		}
		
		else if(i == '\n')
		{
			first_call = false;
			num_cell_refined_ny++;
		}	
	}
	
	fclose(file_U);
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

void write_result_interpolation_to_file (const char filename[], double **A, const double dx, const double dy, const uint32_t nx, const uint32_t ny)
{
	FILE *file = fopen(filename,"w+");	
	
	for (uint32_t i = 0; i < ny; i++)	
	{
		double y = i*dy;
		for (uint32_t j = 0; j < nx; j++)
		{	
			double x = j*dx;
			fprintf(file, "%u %u %g %g %g\n", i, j, x, y, A[i][j]);	
		}
	}
	
	fclose(file);
}

void write_result_to_vtk (const char filename[], double **AA, const double dx, const double dy, const uint32_t nx, const uint32_t ny)
{
	const double dz = 5.0;
	const double CM_TO_UM = 1.0E+04;
	
	FILE *file = fopen(filename,"w+");
	
	fprintf(file,"# vtk DataFile Version 4.2\n");
	fprintf(file,"vtk output\n");
	fprintf(file,"ASCII\n");
	fprintf(file,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(file,"POINTS %u float\n",(nx+1)*(ny+1)*2);
	for (uint32_t k = 0; k < 2; k++)
	{
		double z = k*dz;
		for (uint32_t i = 0; i <= ny; i++)
		{
			double y = i*dy * CM_TO_UM;
			for (uint32_t j = 0; j <= nx; j++)
			{
				double x = j*dx * CM_TO_UM;

				fprintf(file,"%g %g %g\n",x,y,z);

			}
		}	
	}
	fprintf(file,"CELLS %u %u\n",(nx*ny),(nx*ny)*9);
	for (uint32_t i = 0; i < ny; i++)
	{
		for (uint32_t j = 0; j < nx; j++)
		{
			uint32_t k = i*(ny+1) + j;

			uint32_t indexes[8];
			indexes[0] = k;
			indexes[1] = k+1;
			indexes[2] = k+(nx+1)+1;
			indexes[3] = k+(nx+1);
			indexes[4] = k+((nx+1)*(ny+1));
			indexes[5] = k+((nx+1)*(ny+1))+1;
			indexes[6] = k+((nx+1)*(ny+1))+(nx+1)+1;
			indexes[7] = k+((nx+1)*(ny+1))+(nx+1);

			fprintf(file,"8 %u %u %u %u %u %u %u %u\n",indexes[0],indexes[1],indexes[2],indexes[3],indexes[4],indexes[5],indexes[6],indexes[7]);
		}
	}
	fprintf(file,"CELL_TYPES %u\n",(nx*ny));
	for (uint32_t i = 0; i < ny; i++)
	{
		for (uint32_t j = 0; j < nx; j++)
		{
			fprintf(file,"12\n");
		}
	}
	fprintf(file,"CELL_DATA %u\n",(nx*ny));
	fprintf(file,"SCALARS Scalars_ float\n");
	fprintf(file,"LOOKUP_TABLE default\n");
	for (uint32_t i = 0; i < ny; i++)
		for (uint32_t j = 0; j < nx; j++)
			fprintf(file,"%g\n",AA[i][j]);
	
	fclose(file);
}
