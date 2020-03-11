#include "pde_solver.h"

struct solver_pde* new_solver_pde()
{	
	struct solver_pde *pde_data = (struct solver_pde*)malloc(sizeof(struct solver_pde));
	
	return pde_data;
}

void configure_pde_solver(struct solver_pde *s, const double sigma_x, const double sigma_y,\
						const double dt_pde, const double dx_pde, const double dy_pde, const double length_cell_dx, const double length_cell_dy,const uint32_t nx, const uint32_t ny, const uint32_t nt,const char method_name[])
{
	s->dt_pde = dt_pde;
	s->dx_pde = dx_pde;
	s->dy_pde = dy_pde;
	s->sigma_x = sigma_x;
	s->sigma_y = sigma_y;
	s->nx = nx;
	s->ny = ny;
	s->nt = nt;
	s->L_x = length_cell_dx;
	s->L_y = length_cell_dy;
	
	uint32_t size = strlen(method_name);
	s->method_name = (char*)malloc(sizeof(char)*size);
	strcpy(s->method_name,method_name);

	if(strcmp(s->method_name,"DM") == 0 || strcmp(s->method_name,"DDM") == 0)
	{
		s->kappa_x = calculate_kappa(s->L_x,s->dx_pde);
		s->kappa_y = calculate_kappa(s->L_y,s->dy_pde);
	}
	
	else if(strcmp(s->method_name,"CM") == 0)
	{
		
		s->kappa_x = 0.0;
		s->kappa_y = 0.0;
	}
	
	else if(strcmp(s->method_name,"QCM") == 0)
	{
		s->kappa_x = calculate_kappa_QCM(s->L_x,s->dx_pde);
		s->kappa_y = calculate_kappa_QCM(s->L_y,s->dy_pde);			
	}
	

	
}

void allocate_pde_vectors (struct solver_pde *s)
{
	uint32_t nx = s->nx;
	uint32_t ny = s->ny;
	
	s->au_x = (double*)malloc(sizeof(double)*nx);
	s->bu_x = (double*)malloc(sizeof(double)*nx);
	s->cu_x = (double*)malloc(sizeof(double)*nx);
	s->Aux_UX = (double*)malloc(sizeof(double)*nx);
	
	s->av_x = (double*)malloc(sizeof(double)*nx);
	s->bv_x = (double*)malloc(sizeof(double)*nx);
	s->cv_x = (double*)malloc(sizeof(double)*nx);
	s->Aux_VX = (double*)malloc(sizeof(double)*nx);
	
	s->au_y = (double*)malloc(sizeof(double)*ny);
	s->bu_y = (double*)malloc(sizeof(double)*ny);
	s->cu_y = (double*)malloc(sizeof(double)*ny);
	s->Aux_UY = (double*)malloc(sizeof(double)*ny);
	
	s->av_y = (double*)malloc(sizeof(double)*ny);
	s->bv_y = (double*)malloc(sizeof(double)*ny);
	s->cv_y = (double*)malloc(sizeof(double)*ny);
	s->Aux_VY = (double*)malloc(sizeof(double)*ny);
	
	s->rhsu_x = (double*)malloc(sizeof(double)*nx);
	s->rhsu_y = (double*)malloc(sizeof(double)*ny);
	s->rhsv_x = (double*)malloc(sizeof(double)*nx);
	s->rhsv_y = (double*)malloc(sizeof(double)*ny);
		
}

void calculate_parabolic_LHS (struct solver_pde *s)
{
	uint32_t nx = s->nx;
	uint32_t ny = s->ny;
	double dt_pde = s->dt_pde;
	double dx = s->dx_pde;
	double dy = s->dy_pde;
	double sigma_x = s->sigma_x;
	double sigma_y = s->sigma_y;
	
	double L_x = s->L_x;
	double L_y = s->L_y;
	
	double kappa_x = s->kappa_x;
	double kappa_y = s->kappa_y;
	
	// TODO
	// dx_grid	--> ddm
	// L_tot_x	--> ddm
	
	// Set the matrix for the U variable in the x axis
	double *au_x = s->au_x;
	double *bu_x = s->bu_x;
	double *cu_x = s->cu_x;
	set_parabolic_matrix(au_x,bu_x,cu_x,sigma_x,dt_pde,dx,kappa_x,nx);
	
	// Set the matrix for the U variable in the y axis
	double *au_y = s->au_y;
	double *bu_y = s->bu_y;
	double *cu_y = s->cu_y;
	set_parabolic_matrix(au_y,bu_y,cu_y,sigma_y,dt_pde,dy,kappa_y,ny);
	
	// Set the matrix for the V variable in the x axis
	double *av_x = s->av_x;
	double *bv_x = s->bv_x;
	double *cv_x = s->cv_x;
	set_parabolic_matrix(av_x,bv_x,cv_x,sigma_x,dt_pde,dx,kappa_x,nx);
	
	// Set the matrix for the V variable in the y axis
	double *av_y = s->av_y;
	double *bv_y = s->bv_y;
	double *cv_y = s->cv_y;
	set_parabolic_matrix(av_y,bv_y,cv_y,sigma_y,dt_pde,dy,kappa_y,ny);
	
	
}

double calculate_kappa(double L_h, double h)
{
	
	double result;
	
	result = (L_h*L_h*L_h*L_h - h*h*h*h)/(12.0*(L_h*L_h)*(h*h));
	
	return result;
	
}

double calculate_kappa_QCM(double L_h, double h)
{
	
	double result;
	
	result = (L_h*L_h)/(12.0*(h*h));
	
	return result;
	
}

void set_parabolic_matrix (double *a, double *b, double *c, const double sigma, const double dt, const double h, const double kappa, const uint32_t n)
{
	/*
		a: V[i-1] coefficient - inferior diagonal
		b: V[i] coefficient - main diagonal
		c: V[i+1] coefficient - superior diagonal
	*/	
	
	for (uint32_t i = 0; i < n; i++)
	{
		// Left border
		if (i == 0)
		{
			a[i] = 0.0;
			//~ b[i] = (dt/(2.0*h*h))*(sigma) + 1.0 + calculate_kappa(L_h,h);
			//~ c[i] = -(dt/(2.0*h*h))*(sigma) - calculate_kappa(L_h,h);
			//b[i] = (dt/2.0)*(transm(i,i+1,sigma_0,dx)) + 1.0 + (transm_k(i,i+1,dx,L));	// --> ddm
			//c[i] = -((dtime/2.0)*transm(i,i+1,sigma_0,dx) + transm_k(i,i+1,dx,L));		// --> ddm
			b[i] = (dt/(2.0*h*h))*(sigma) + 1.0 + kappa;
			c[i] = -(dt/(2.0*h*h))*(sigma) - kappa;			
		}
		// Right border
		else if (i == (n-1))
		{
			//a[i] = -((dtime/2.0)*transm(i,i-1,sigma_0,dx) + transm_k(i,i-1,dx,L));		// --> ddm
			//b[i] =  (dtime/2.0)*(transm(i,i-1,sigma_0,dx)) + 1.0 + (transm_k(i,i-1,dx,L));// --> ddm
			//~ a[i] = -(dt/(2.0*h*h))*(sigma) - calculate_kappa(L_h,h);
			//~ b[i] = (dt/(2.0*h*h))*(sigma) + 1.0 + calculate_kappa(L_h,h);
			a[i] = -(dt/(2.0*h*h))*(sigma) - kappa;
			b[i] = (dt/(2.0*h*h))*(sigma) + 1.0 + kappa;
			c[i] = 0.0;
		}
		else
		{
			//a[i] = -((dtime/2.0)*transm(i,i-1,sigma_0,dx) + transm_k(i,i-1,dx,L));
			//b[i] =  (dtime/2.0)*(transm(i,i+1,sigma_0,dx) + transm(i,i-1,sigma_0,dx)) + 1.0 + (transm_k(i,i+1,dx,L)+transm_k(i,i-1,dx,L));
			//c[i] = -((dtime/2.0)*transm(i,i+1,sigma_0,dx) + transm_k(i,i+1,dx,L));
			//~ a[i] = -(dt/(2.0*h*h))*(sigma) - calculate_kappa(L_h,h);
			//~ b[i] =  (dt/(2.0*h*h))*( (sigma) + (sigma) ) + 1.0 + 2.0*calculate_kappa(L_h,h);
			//~ c[i] = -(dt/(2.0*h*h))*(sigma) - calculate_kappa(L_h,h);
			a[i] = -(dt/(2.0*h*h))*(sigma) - kappa;
			b[i] =  (dt/(2.0*h*h))*( (sigma) + (sigma) ) + 1.0 + 2.0*kappa;
			c[i] = -(dt/(2.0*h*h))*(sigma) - kappa;
		}
	}
}

void solve_pde (struct solver_pde *the_pde_solver, double **U, double **V)
{
	uint32_t nx = the_pde_solver->nx;
	uint32_t ny = the_pde_solver->ny;

	double kappa_x = the_pde_solver->kappa_x;
	double kappa_y = the_pde_solver->kappa_y;

	double L_x = the_pde_solver->L_x;
	double L_y = the_pde_solver->L_y;
	
	double dx = the_pde_solver->dx_pde;
	double dy = the_pde_solver->dy_pde;
	
	// TODO: Calculate kappa
	// double kappa_x, kappa_y;
	// if (method == "DM" || method == "DDM")
	// 		kappa_x = calculate_kappa(L_x,dx)
	// 		kappa_y = calculate_kappa(L_y,dy)
	// else if (method == "CM"
	// 		kappa_x = 0.0
	//      kappa_y = 0.0
	// else if (method == "QSM")
	//		kappa_x = calcula_kappa_qsm(L_x,dx)
	//		kappa_y = calcula_kappa_qsm(L_y,dy)
	
	// Activator U --> X axis
	double *au_x = the_pde_solver->au_x;
	double *bu_x = the_pde_solver->bu_x;
	double *cu_x = the_pde_solver->cu_x;
	double *rhsu_x = the_pde_solver->rhsu_x;
	solve_parabolic(au_x,bu_x,cu_x,rhsu_x,U,nx,ny,'x', kappa_x);
	
	// Inactivator V --> X axis
	double *av_x = the_pde_solver->av_x;
	double *bv_x = the_pde_solver->bv_x;
	double *cv_x = the_pde_solver->cv_x;
	double *rhsv_x = the_pde_solver->rhsv_x;
	solve_parabolic(av_x,bv_x,cv_x,rhsv_x,V,nx,ny,'x', kappa_x);

	// Activator U --> Y axis
	double *au_y = the_pde_solver->au_y;
	double *bu_y = the_pde_solver->bu_y;
	double *cu_y = the_pde_solver->cu_y;
	double *rhsu_y = the_pde_solver->rhsu_y;
	solve_parabolic(au_y,bu_y,cu_y,rhsu_y,U,nx,ny,'y', kappa_y);
	
	// Inactivator V --> Y axis
	double *av_y = the_pde_solver->av_y;
	double *bv_y = the_pde_solver->bv_y;
	double *cv_y = the_pde_solver->cv_y;
	double *rhsv_y = the_pde_solver->rhsv_y;
	solve_parabolic(av_y,bv_y,cv_y,rhsv_y,V,nx,ny,'y', kappa_y);
}

void solve_parabolic (const double *a, const double *b, const double *c, double *rhs, double **A, const uint32_t nx, const uint32_t ny, const char axis, const double kappa)
{
	// TODO: If it is a refined mesh we need to allocate here ...
	double aux_x[nx], aux_y[ny];
	
	if (axis == 'x')
	{
		
		for (uint32_t i = 0; i < ny; i++)
		{
			for (uint32_t j = 0; j < nx; j++)
			{
				// Left side
				if (j == 0)
				{
					//rhsX[j] = -(transm_k(j,j+1,dx_grid,L_tot_x))*U[i][j+1] + (1.0 + transm_k(j,j+1,dx_grid,L_tot_x))*U[i][j]; // ---> ddm
					//~ rhs[j] =  +(1.0)*A[i][j];
					//~ rhs[j] = -(calculate_kappa(L_h,h))*A[i][j+1] +(1.0 + calculate_kappa(L_h,h))*A[i][j];
					//~ rhs[j] = -(calculate_kappa(L_h,h))*A[i][j+1] +(1.0 + calculate_kappa(L_h,h))*A[i][j];
					rhs[j] = -(kappa)*A[i][j+1] +(1.0 + kappa)*A[i][j];
				}
				// Right side
				else if (j == (nx-1))
				{
					//rhsX[j] = -(transm_k(j,j-1,dx_grid,L_tot_x))*U[i][j-1] + (1.0 + transm_k(j,j-1,dx_grid,L_tot_x))*U[i][j]; // --> ddm
					//~ rhs[j] = +(1.0)*A[i][j];
					//~ rhs[j] = -(calculate_kappa(L_h,h))*A[i][j-1] +(1.0 + calculate_kappa(L_h,h))*A[i][j];
					rhs[j] = -(kappa)*A[i][j-1] +(1.0 + kappa)*A[i][j];
				}
				else
				{
					//rhsX[j] = -(transm_k(j,j+1,dx_grid,L_tot_x))*U[i][j+1] + (1.0 + (transm_k(j,j+1,dx_grid,L_tot_x) + transm_k(j,j-1,dx_grid,L_tot_x)))*U[i][j] - transm_k(j,j-1,dx_grid,L_tot_x)*U[i][j-1]; // --> ddm
					//~ rhs[j] =  +(1.0)*A[i][j];
					//~ rhs[j] = -(calculate_kappa(L_h,h))*A[i][j+1] +(1.0 + 2.0*(calculate_kappa(L_h,h)))*A[i][j] -(calculate_kappa(L_h,h))*A[i][j-1];
					rhs[j] = -(kappa)*A[i][j+1] +(1.0 + 2.0*(kappa))*A[i][j] -(kappa)*A[i][j-1];
				}
				aux_x[j] = rhs[j];
			}
			
			solve_tridiagonal_system(a,b,c,rhs,aux_x,nx);
			
			for (uint32_t j = 0 ; j < nx; j++)
			{			
				A[i][j] = aux_x[j];
			}
		}
	}
	else if (axis == 'y')
	{
		for (uint32_t i = 0; i < nx; i++)
		{
			for (uint32_t j = 0; j < ny; j++)
			{
				// Left side
				if (j == 0)
				{
					//~ rhs[j] = +(1.0)*A[j][i];
					//~ rhs[j] = -(calculate_kappa(L_h,h))*A[j+1][i] +(1.0 + calculate_kappa(L_h,h))*A[j][i];
					rhs[j] = -(kappa)*A[j+1][i] +(1.0 + kappa)*A[j][i];
				}
				// Right side
				else if (j == (ny-1))
				{
					//~ rhs[j] = +(1.0)*A[j][i];
					//~ rhs[j] = -(calculate_kappa(L_h,h))*A[j-1][i]+(1.0 + calculate_kappa(L_h,h))*A[j][i];
					rhs[j] = -(kappa)*A[j-1][i]+(1.0 + kappa)*A[j][i];
				}
				else
				{
					//~ rhs[j] = +(1.0)*A[j][i];
					//~ rhs[j] = -(calculate_kappa(L_h,h))*A[j+1][i] + (1.0+ 2.0*calculate_kappa(L_h,h))*A[j][i] -(calculate_kappa(L_h,h))*A[j-1][i];
					rhs[j] = -(kappa)*A[j+1][i] + (1.0+ 2.0*kappa)*A[j][i] -(kappa)*A[j-1][i];
				}
				
				aux_y[j] = rhs[j];
			}
			
			solve_tridiagonal_system(a,b,c,rhs,aux_y,ny);
			
			for (uint32_t j = 0 ; j < ny; j++)
			{			
				A[j][i] = aux_y[j];
			}
		}
	}
	else
	{
		fprintf(stderr,"[-] ERROR! Invalid axis for 'solve_parabolic'\n");
		exit(EXIT_FAILURE);
	}
	
}

void solve_tridiagonal_system (const double *a, const double *b, const double *c, double *rhs, double *x, const uint32_t n)
{
	double *aux_c, *aux_rhs;
	aux_c = (double*)malloc(sizeof(double)*n);
	aux_rhs = (double*)malloc(sizeof(double)*n);
	
	if(b[0] == 0) 
	{
		fprintf(stderr,"[-] ERROR! The matrix is singular\n");
		exit(EXIT_FAILURE);
	}
	else 
	{
		aux_c[0] = c[0]/b[0];
		aux_rhs[0] = rhs[0]/b[0];
	}
	
	for (uint32_t i = 1; i < n; i++) 
	{
		if( (b[i] - aux_c[i-1] * a[i]) == 0) 
		{
			fprintf(stderr,"[-] ERROR! The matrix is singular\n");
			exit(EXIT_FAILURE);
		}
		else 
		{
			double id = 1.0/(b[i] - aux_c[i-1] * a[i]);
			aux_c[i] = c[i]*id; 
			aux_rhs[i] = (rhs[i] - aux_rhs[i-1] * a[i]) * id;
		}
	}
	
	/* Retro-substitition */
	x[n-1] = aux_rhs[n-1];
	for (int i = n-2; i >= 0; i--)
	{
		x[i] = aux_rhs[i] - aux_c[i] * x[i+1];
	}
	
	free(aux_c);
	free(aux_rhs);
}
