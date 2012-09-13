
//ebss:
#include<common/parameters/HamiltonianParameters.hpp>
#include<common/common.hpp>

//petsc:
#include<petsc.h>

//stl:
#include<vector>
#include<sstream>
#include<string>
#include<iostream>


int main(int argc, char ** argv)
{
    PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

    HamiltonianParameters<PetscReal> *params = new HamiltonianParameters<PetscReal>(MPI_COMM_WORLD);
    params->print_parameters();

    delete params;
    PetscFinalize();
    return 0;
}
