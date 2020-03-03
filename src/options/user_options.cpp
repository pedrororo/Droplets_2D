#include "user_options.h"

struct user_options* new_user_options (int argc, char *argv[])
{
    struct user_options *result = (struct user_options*)malloc(sizeof(struct user_options));

	result->spiral_time = -1;
	result->save_state_filename = NULL;
	result->input_filename_U = NULL;
	result->input_filename_V = NULL;

    read_config_file(result,argv[1]);
    
    return result;
}

void free_user_options (struct user_options *options)
{   
    free(options);
}

void read_config_file (struct user_options *options, const char filename[])
{
    printf("%s\n",PRINT_LINE);
    printf("[user_options] Reading configuration file:> \"%s\"\n",filename);

    // Open the config file for reading
    FILE *file = fopen(filename,"r");
    if (!file)
    {
        fprintf(stderr,"[-] Error reading configuration file '%s'\n",filename);
        exit(EXIT_FAILURE);
    }
    
    // Here we parse the config file
    if(ini_parse(filename, parse_config_file, options) < 0) 
    {
        fprintf(stderr, "Error: Can't load the config file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    printf("%s\n",PRINT_LINE);

    fclose(file);

	// DEBUG
    print_user_options(options);
}

int parse_config_file(void *user, const char *section, const char *name, const char *value) 
{
    struct user_options *pconfig = (struct user_options *)user;

    if (SECTION_STARTS_WITH(MAIN_SECTION))
    {
        if (MATCH_NAME("dt_pde"))
        {
            pconfig->dt_pde = strtof(value, NULL);
        }
        else if (MATCH_NAME("simulation_time"))
        {
            pconfig->simulation_time = strtof(value, NULL);
        }
        else if (MATCH_NAME("spiral_time"))
        {
            pconfig->spiral_time = strtof(value, NULL);
        }
    }
    else if (SECTION_STARTS_WITH(SAVE_STATE_SECTION))
    {
        if (MATCH_NAME("save_rate"))
        {
            pconfig->save_rate = strtol(value,NULL,10);
        }
        else if (MATCH_NAME("save_state_filename"))
        {
            pconfig->save_state_filename = strdup(value);
        }
    }
    else if (SECTION_STARTS_WITH(LOAD_STATE_SECTION))
    {
        if (MATCH_NAME("input_filename_U"))
        {
            pconfig->input_filename_U = strdup(value);
        }
        else if (MATCH_NAME("input_filename_V"))
        {
            pconfig->input_filename_V = strdup(value);
        }
    }
    else if (SECTION_STARTS_WITH(UPDATE_CHEMICAL_SECTION))
    {
        if (MATCH_NAME("main_function"))
        {
            pconfig->update_chemical_function_name = strdup(value);
        }
    }
    else if (SECTION_STARTS_WITH(SAVE_RESULT_SECTION))
    {

        if (MATCH_NAME("output_dir"))
        {
            pconfig->output_dir = strdup(value);
        }   
        else if (MATCH_NAME("print_rate"))
        {
            pconfig->print_rate = strtol(value,NULL,10);
        }
    }
    else if (SECTION_STARTS_WITH(ASSEMBLY_MATRIX_SECTION))
    {

        if (MATCH_NAME("sigma_x"))
        {
            pconfig->sigma_x = strtof(value, NULL);
        }
        else if (MATCH_NAME("sigma_y"))
        {
            pconfig->sigma_y = strtof(value, NULL);
        }
        else if (MATCH_NAME("sigma_z"))
        {
            pconfig->sigma_z = strtof(value, NULL);
        }   
        else if (MATCH_NAME("main_function"))
        {
            pconfig->assembly_matrix_function_name = strdup(value);
        }
    }
    else if (SECTION_STARTS_WITH(LINEAR_SYSTEM_SECTION))
    {
        if (MATCH_NAME("main_function"))
        {
            pconfig->linear_system_function_name = strdup(value);
        }
    }
    else if (SECTION_STARTS_WITH(DOMAIN_SECTION))
    {

        if (MATCH_NAME("start_dx"))
        {
            pconfig->start_dx = strtof(value, NULL);
        }
        else if (MATCH_NAME("start_dy"))
        {
            pconfig->start_dy = strtof(value, NULL);
        }
        else if (MATCH_NAME("start_dz"))
        {
            pconfig->start_dz = strtof(value, NULL);
        }   
        else if (MATCH_NAME("side_length"))
        {
            pconfig->side_length = strtof(value, NULL);
        }
        else if (MATCH_NAME("main_function"))
        {
            pconfig->domain_function_name = strdup(value);
        }
    }
    else if (SECTION_STARTS_WITH(ODE_SOLVER_SECTION))
    {

        if (MATCH_NAME("dt_ode"))
        {
            pconfig->dt_ode = strtof(value, NULL);
        }
    }
    else if (SECTION_STARTS_WITH(STIM_S1_SECTION))
    {
		// TODO
    }
    else if (SECTION_STARTS_WITH(STIM_S2_SECTION))
    {
        // TODO
    }

    return 1;
}

void print_user_options (struct user_options *options)
{
    printf("********************* user_options *********************\n");
    printf("dt_pde = %g\n",options->dt_pde);
    printf("simulation_time = %g\n",options->simulation_time);
    printf("--------------------------------------------------------\n");
    printf("save_rate = %u\n",options->save_rate);
    printf("save_state_filename = %s\n",options->save_state_filename);
    printf("--------------------------------------------------------\n");
    printf("update_chemical_function_name = %s\n",options->update_chemical_function_name);
    printf("--------------------------------------------------------\n");
    printf("print_rate = %u\n",options->print_rate);
    printf("output_dir = %s\n",options->output_dir);
    printf("--------------------------------------------------------\n");
    printf("sigma_x = %g\n",options->sigma_x);
    printf("sigma_y = %g\n",options->sigma_y);
    printf("sigma_z = %g\n",options->sigma_z);
    printf("assembly_matrix_function_name = %s\n",options->assembly_matrix_function_name);
    printf("--------------------------------------------------------\n");
    printf("linear_system_function_name = %s\n",options->linear_system_function_name);
    printf("--------------------------------------------------------\n");
    printf("start_dx = %g\n",options->start_dx);
    printf("start_dy = %g\n",options->start_dy);
    printf("start_dz = %g\n",options->start_dz);
    printf("side_length = %lf\n",options->side_length);
    printf("domain_function_name = %s\n",options->domain_function_name);
    printf("--------------------------------------------------------\n");
    printf("dt_ode = %lf\n",options->dt_ode);
    printf("--------------------------------------------------------\n");
    // TODO (Stimuli ...)
    printf("********************************************************\n");
}
