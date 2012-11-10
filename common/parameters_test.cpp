#include<common/new_parameters/BasisParameters.hpp>
#include<common/new_parameters/HamiltonianParameters.hpp>
//#include<common/common.hpp>
#include<vector>
#include<petsc.h>

int
main (int argc, const char* argv[])
{
    int ac = argc;
    char** av = new char*[argc];
    for (size_t i = 0; i < argc; i++)
    {
        av[i] = new char[strlen(argv[i])];
        std::copy(argv[i], argv[i] + strlen(argv[i]), av[i]);
    }

    PetscInitialize(&ac,&av,PETSC_NULL,PETSC_NULL);

    BasisParameters<double, double> p(argc, argv, PETSC_COMM_WORLD);
    HamiltonianParameters h(argc, argv , PETSC_COMM_WORLD);

    std::cout << p.print();
}
