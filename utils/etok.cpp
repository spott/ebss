//ebss:
#include<common/common.hpp>
#include<common/parameters/KParameters.hpp>
#include<common/coulomb/complex_functions.H>

//petsc:
#include<petsc.h>

int main(int argc, const char ** argv)
{
    int ac = argc;
    char** av = new char*[argc];
    for (size_t i = 0; i < argc; i++)
    {
        av[i] = new char[strlen(argv[i])+1];
        std::copy(argv[i], argv[i] + strlen(argv[i])+1, av[i]);
    }
    PetscInitialize(&ac, &av, PETSC_NULL, PETSC_NULL);

    KParameters *kparams = new KParameters(argc, argv, MPI_COMM_WORLD);
    //BasisParameters *bparams = new BasisParameters(argc, argv, MPI_COMM_WORLD);
    if (kparams->rank() == 0)
    {
        std::cout << kparams->print();
        //std::cout << bparams->print();
    }


