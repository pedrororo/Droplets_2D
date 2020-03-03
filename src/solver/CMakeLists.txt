set(SOURCE_FILES_SOLVER pde_solver.h pde_solver.cpp
						ode_solver.h ode_solver.cpp
						droplets_solver.h droplets_solver.cpp)

ADD_LIBRARY(solver STATIC ${SOURCE_FILES_SOLVER})

if(NOT MSVC)
    set_property(TARGET solver PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
