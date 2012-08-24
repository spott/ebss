#include<iostream>

#include<common/parameters.hpp>
#include<petsc.h>

int main(int argc, char *argv[])
{
   PetscErrorCode ierr;

   PetscInitialize(&argc,&argv);
   MPI_Comm world = MPI_COMM_WORLD;
   parameters params(world, BASIS);
   params.print_parameters();
   std::cout << "a" << std::endl;
   params.save_parameters("stuff");
   std::cout << "b" << std::endl;

   parameters p2(world, BASIS);
   p2.init_from_file("stuff");
   std::cout << "c" << std::endl;
   ierr = p2.print_parameters(); CHKERRQ(ierr);
   std::cout << p2.basis_size() << std::endl;

   PetscFinalize();
   return 0;
}
