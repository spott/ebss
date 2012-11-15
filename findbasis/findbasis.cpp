#include<common/parameters/BasisParameters.hpp>
#include<findbasis/numerov.hpp>
#include<findbasis/finite_difference.hpp>
//#include<common/special/bspline.hpp>

#include<slepc.h>

struct sae_param {
    int p;
    double c;
    double beta;
};

template <typename T>
T hydrogen_pot(T r)
{
   return -1./r;
}

template <typename T>
T neon_pot(T r, T Z, T N, std::vector<sae_param> atom)
{
    T a = 0;
    for( auto p: atom )
        a += p.c * std::pow(r, p.p) * std::exp(- p.beta * r);
    a *= (N-1);
    a += 1 - N + Z;
    a *= -1/r;
   return a;
}


typedef double scalar;

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


    std::vector< sae_param > neon{
        {0, .46087879, 4.68014471},
        {0, 0.53912121, 2.41322960},
        {1, 0.42068967, 5.80903874},
        {1, 0.47271993, 2.90207510},
        {2, -1.12569309, 4.51696279},
        {2, 1.29942636, 3.06518063}};


    //for (size_t i = 0; i < 10000; i++)


    //call function to find all the energy states here:
    //numerov::find_basis_set<scalar>( [neon](scalar r) {return neon_pot<scalar>(r, 10, 10, neon);}, params);
    
    finite_difference::find_basis<2, scalar>( [neon](scalar r) {return neon_pot<scalar>(r, 10,10, neon);}, 
                                              params);
    //finite_difference::find_basis<2, scalar>( hydrogen_pot<scalar>, 
                                              //params);

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
