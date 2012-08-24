#include<common/parameters.hpp>

#include<petsc.h>

int
main (int argc, char** argv)
{
    PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
    BasisParameters p(PETSC_COMM_WORLD);

}
