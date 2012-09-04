#include<common/parameters.hpp>
#include<common/common.hpp>
#include<vector>
#include<petsc.h>

int
main (int argc, char** argv)
{
    PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
    BasisParameters p(PETSC_COMM_WORLD);

    //std::vector<double> grid = common::import_vector_binary<double>(p.grid_file());

    std::vector<PetscScalar> *grid = p.grid();
    for (size_t i = 0; i < grid->size(); i++)
        std::cout << grid->at(i) << ", ";
    std::cout << std::endl;
    p.print_parameters();
}
