
//ebss:
#include<common/parameters/HamiltonianParameters.hpp>
#include<common/common.hpp>
#include<common/math.hpp>

//petsc:
#include<petsc.h>

//stl:
#include<vector>
#include<sstream>
#include<string>
#include<iostream>


int main(int argc, const char ** argv)
{
    int ac = argc;
    char** av = new char*[argc];
    for (size_t i = 0; i < argc; i++)
    {
        av[i] = new char[strlen(argv[i])];
        std::copy(argv[i], argv[i] + strlen(argv[i]), av[i]);
    }
    PetscInitialize(&ac, &av, PETSC_NULL, PETSC_NULL);

    HamiltonianParameters<PetscReal> *params = new HamiltonianParameters<PetscReal>(argc, argv, MPI_COMM_WORLD);
    if (params->rank() == 0)
        std::cout << params->print();

    std::vector<PetscReal> *grid = params->basis_parameters()->grid();
    auto prototype = params->prototype();

    std::function<bool (int, int)> dipole_selection_rules = [prototype](int i, int j) {
        return (
                (std::abs(prototype[i].l - prototype[j].l) == 1) || i == j
               );
    };

    std::function<PetscScalar (int, int)> findvalue = [prototype,params,grid](int i, int j)->PetscScalar{
        if (i == j)
            return 0.0;
        PetscScalar radial = math::integrateGrid(
                common::import_vector_binary<PetscReal>(params->basis_parameters()->basis_function_filename(prototype[i])),
                common::import_vector_binary<PetscReal>(params->basis_parameters()->basis_function_filename(prototype[j])),
                *grid);
        PetscScalar angular = math::CGCoefficient<PetscScalar>(prototype[i],prototype[j]);
       return radial * angular;
    };

   Mat H = common::populate_matrix(*params, 
                                    dipole_selection_rules, 
                                    findvalue, 
                                    prototype.size(),
                                    prototype.size());


   params->write_dipole_matrix(H);
   params->write_energy_eigenvalues();

   params->save_parameters();
   delete params;
   PetscFinalize();
   return 0;
}
