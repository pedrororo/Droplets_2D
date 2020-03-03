//
// Created by bergolho on 12/02/19.
//

#ifndef USER_OPTIONS_H
#define USER_OPTIONS_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>

#include "../utils/utils.h"
#include "../ini_parser/ini.h"
#include "../ini_parser/ini_file_sections.h"

//#include "cost_function_config.h"
//#include "local_optimization_config.h"

static const uint32_t MAX_FILENAME_SIZE = 200;

struct user_options
{
    
    double dt_pde;
    double simulation_time;
    
    uint32_t save_rate;
    char *save_state_filename;
    
    char *input_filename_U;
    char *input_filename_V;

	char *update_chemical_function_name;
	
	uint32_t print_rate;
	char *output_dir;
	
	double sigma_x;
	double sigma_y;
	double sigma_z;
	char *assembly_matrix_function_name;
	
	char *linear_system_function_name;
	
	double start_dx;
	double start_dy;
	double start_dz;
	double side_length;
	char *domain_function_name;
	
	double dt_ode;
	
	double spiral_time;
	
	//struct stimuli_config *stim_config_s1;
	//struct stimuli_config *stim_config_s2;

};

struct user_options* new_user_options (int argc, char *argv[]);
void free_user_options (struct user_options *options);

void read_config_file (struct user_options *options, const char filename[]);

int parse_config_file(void *user, const char *section, const char *name, const char *value);

void print_user_options (struct user_options *options);

#endif
