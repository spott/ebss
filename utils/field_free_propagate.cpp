//ebss:
#include<common/parameters/HamiltonianParameters.hpp>
#include<common/math.hpp>

//stl:
#include<iostream>


int
main(const int argc, const char ** argv)
{
    int ac = argc;
    char** av = new char*[argc];
    for (size_t i = 0; i < argc; i++)
    {
        av[i] = new char[strlen(argv[i])+1];
        std::copy(argv[i], argv[i] + strlen(argv[i])+1, av[i]);
    }
    PetscInitialize(&ac, &av, PETSC_NULL, PETSC_NULL);

    PetscBool flg = PETSC_FALSE;
    char bagname[PETSC_MAX_PATH_LEN];
    PetscOptionsGetString(PETSC_NULL, "-hamiltonian_config", bagname, PETSC_MAX_PATH_LEN, &flg);

    if (!flg)
    {
        std::cerr << "I need a hamiltonian to propagate. (-hamiltonian_config )" << std::endl;
        PetscFinalize();
        return 0;
    }


    //get dt and max value for time prop:
    PetscReal dt;
    PetscReal max_t;
    PetscOptionsGetReal(PETSC_NULL, <#"-dt", &dt, &flg);
    PetscOptionsGetReal(PETSC_NULL, <#"-max_t", &maxt_t, &flg);
    //PetscOptionsGetReal(PETSC_NULL, "-dt", bagname, PETSC_MAX_PATH_LEN, &flg);
    //if (!flg)
    //{
        //std::cerr << "I need a hamiltonian to propagate. (-hamiltonian_config )" << std::endl;
        //PetscFinalize();
        //return 0;
    //}

    //HamiltonianParameters<PetscReal> *params = new HamiltonianParameters<PetscReal>(MPI_COMM_WORLD, std::string(bagname) );

    //return 0;
}
