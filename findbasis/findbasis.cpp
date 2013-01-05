#include<common/parameters/BasisParameters.hpp>
#include<findbasis/numerov.hpp>
#include<common/math.hpp>

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
    if (params->rank() == 0) 
        std::cout << params->print();

    sae<scalar> hydrogen = {potassium_params, 1, 1, -.5};


    //call function to find all the energy states here:
    numerov::find_basis_set<scalar>( [](scalar r) {return pot<scalar>(r);}, params);

    //write out parameters:
    params->save_parameters();

    //delete the params...
    delete params;
    SlepcFinalize();
    return 0;
}
