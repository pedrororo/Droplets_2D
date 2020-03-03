set(SOURCE_FILES_PDE pde_solver.h pde_solver.cpp)

ADD_LIBRARY(pde_solver STATIC ${SOURCE_FILES_PDE})

if(NOT MSVC)
    set_property(TARGET pde_solver PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
