#include<common/parameters/BasisParameters.hpp>
#include<findbasis/numerov.hpp>
//#include<common/special/bspline.hpp>

#include<slepc.h>

template <typename T>
T pot(T r)
{
   return -1./r;
}

typedef long double scalar;

int main(int argc, char **argv)
{
   SlepcInitialize(&argc, &argv, PETSC_NULL, "FindBasis - Find a numerical basis");
   MPI_Comm world = MPI_COMM_WORLD;

   BasisParameters<scalar> *params = new BasisParameters<scalar>(world);

   // print out the parameters
   params->print_parameters();

   //call function to find all the energy states here:
   numerov::find_basis_set<scalar>(&pot, params);

   //delete the params... this is important! it writes out everything.
   delete params;
   SlepcFinalize();
   return 0;
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
 * find_basis_set< ScalarType >( potential function, l_value , BasisParameters)
 *
 * This should be all that is needed.  So using different
 * find_basis functions, SHOULD make this fairly uniform.
 *
 * Different grids should be usable, but the grid should be in the
 * parameters class anyways.
 *
 **************************************/
