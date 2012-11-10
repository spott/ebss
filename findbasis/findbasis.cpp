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

int main(int argc, const char **argv)
{
    int ac = argc;
    char** av = new char*[argc];
    for (size_t i = 0; i < argc; i++)
    {
        av[i] = new char[strlen(argv[i])];
        std::copy(argv[i], argv[i] + strlen(argv[i]), av[i]);
    }

    SlepcInitialize(&ac, &av, PETSC_NULL, "FindBasis - Find a numerical basis");

    BasisParameters<scalar> *params = new BasisParameters<scalar>(argc, argv, PETSC_COMM_WORLD);

    // print out the parameters
    std::cout << params->print();

    //call function to find all the energy states here:
    numerov::find_basis_set<scalar>(&pot, params);

    //write out parameters:
    params->save_parameters();

    //delete the params...
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
