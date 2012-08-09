#include<common/parameters.hpp>

#include<boost/bind.hpp>

#include<findbasis/numerov.hpp>
//#include<common/special/bspline.hpp>

#include<slepc.h>

PetscReal pot(PetscReal r, int l);

int main(int argc, char *argv[])
{

   //PetscInitialize(&argc, &argv);
   SlepcInitialize(&argc, &argv, (char*)0, "");
   MPI_Comm world = MPI_COMM_WORLD;

   parameters params(world, BASIS);
   params.print_parameters();

   //call specific function here:
   numerov::find_basis<PetscReal>(
      params,
      boost::bind<PetscReal>(pot, _1, (int)0));

   PetscFinalize();
   return 0;
}

PetscReal pot(PetscReal r, int l)
{
   return -1/r - l*(l+1)/(r*r);
}

/**************************************
 *
 * We need to design an interface for the 
 * basis finder.  The potential function,
 * sans the angular term, seems most appropriate.
 * It seems unlikely that I will ever be using this code 
 * (without major modification) on non-spherically
 * symmetric code.
 *
 * find_basis< ScalarType >( potential function, l_value , parameters)
 *
 * This should be all that is needed.  So using different
 * find_basis functions, SHOULD make this fairly uniform.
 *
 * Different grids should be usable, but the grid should be in the
 * parameters class anyways.
 *
 **************************************/
