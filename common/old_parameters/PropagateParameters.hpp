#pragma once

//ebss:
#include<common/common.hpp>
#include<common/parameters/Parameters.hpp>

//stl:
#include<sstream>
#include<string>

//petsc:
#include<petsc.h>

class PropagateParameters: public Parameters
{
public:
    typedef struct {
        char hamiltonian_bag_filename[PETSC_MAX_PATH_LEN];
    }
}
