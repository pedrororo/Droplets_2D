//
// Created by sachetto on 12/10/17.
//

#ifndef MONOALG3D_INI_FILE_HEADERS_H
#define MONOALG3D_INI_FILE_HEADERS_H

#define MAIN_SECTION "main"
#define SAVE_STATE_SECTION "save_state"
#define LOAD_STATE_SECTION "load_state"
#define UPDATE_CHEMICAL_SECTION "update_chemical"
#define SAVE_RESULT_SECTION "save_result"
#define ASSEMBLY_MATRIX_SECTION "assembly_matrix"
#define LINEAR_SYSTEM_SECTION "linear_system_solver"
#define DOMAIN_SECTION "domain"
#define ODE_SOLVER_SECTION "ode_solver"
#define STIM_S1_SECTION "stim_s1"
#define STIM_S2_SECTION "stim_s2"

#define MATCH_SECTION_AND_NAME(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
#define MATCH_SECTION(s) strcmp(section, s) == 0
#define MATCH_NAME(v) strcmp(name, v) == 0
#define SECTION_STARTS_WITH(s) strncmp(section, s, strlen(s)) == 0
#define A_STARTS_WITH_B(a ,b) strncmp(a, b, strlen(b)) == 0


#endif //MONOALG3D_INI_FILE_HEADERS_H
