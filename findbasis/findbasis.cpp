#include <common/parameters/BasisParameters.hpp>
#include <findbasis/numerov.hpp>
#include <findbasis/single_active_electron.hpp>
//#include<findbasis/finite_difference.hpp>
//#include<common/special/bspline.hpp>
//#include<common/math.hpp>
#include <functional>

#include <slepc.h>


typedef long double scalar;
int main( int argc, const char** argv )
{
    int    ac = argc;
    char** av = new char*[argc + 1];
    for ( size_t i = 0; i < argc; i++ ) {
        av[i] = new char[strlen( argv[i] )];
        std::copy( argv[i], argv[i] + strlen( argv[i] ), av[i] );
    }
    av[argc] = NULL;

    SlepcInitialize( &ac, &av, PETSC_NULL,
                     "FindBasis - Find a numerical basis" );

    BasisParameters<scalar> params( argc, argv, PETSC_COMM_WORLD );

    // print out the parameters
    if ( params.rank() == 0 ) {
        std::ofstream out( "git_commit.txt" );
        out << GIT_COMMIT << std::endl;
        out.close();
        std::cout << "#Git commit: " << GIT_COMMIT << std::endl;
        std::cout << params.print();
    }

    std::vector<sae_param<scalar>> helium_params{
        {0, 0.33740902, 6.28379042},  {0, 0.66259098, 0.23817624},
        {1, 1.67240025, 3.60759648},  {1, -1.15584925, 1.39690934},
        {2, -1.10884184, 3.64179068}, {2, -0.33200330, 0.81080622},
    };
    sae<scalar> helium = {helium_params, 2, 2, -0.91793786};

    std::vector<sae_param<scalar>> neon_params{
        {0, .46087879, 4.68014471},   {0, 0.53912121, 2.41322960},
        {1, 0.42068967, 5.80903874},  {1, 0.47271993, 2.90207510},
        {2, -1.12569309, 4.51696279}, {2, 1.29942636, 3.06518063}};
    sae<scalar> neon = {neon_params, 10, 10, -30.87114223};

    std::vector<sae_param<scalar>> argon_params{
        {0, -0.14865117, 1.69045518},   {0, 1.14865117, 2.48113492},
        {1, -2.01010923, 83.81780982},  {1, -0.01998004, 0.57781500},
        {2, -20.01846287, 14.79181994}, {2, 0.53596000, 1.83863745}};
    sae<scalar> argon = {argon_params, 18, 18, -114.45726722};

    std::vector<sae_param<scalar>> rubidium_params{
        {0, 0.81691787, 7.83077875}, {0, 0.18308213, 2.75163799},
        {1, 2.53670563, 4.30010258}, {2, -19.56508990, 43.31975597},
        {3, 1.06320272, 2.93818679}, {4, -0.99934358, 4.97097146}};
    sae<scalar> rubidium = {rubidium_params, 37, 37, -541.83186666};

    std::vector<sae_param<scalar>> potassium_params{
        {0, 0.77421694, 11.32240776}, {0, 0.22578306, 2.32391666},
        {1, 6.07532608, 5.49215770},  {2, -17.36457454, 11.17407777},
        {3, 1.95999037, 2.84608178},  {4, -0.13294584, 2.36656543}};
    sae<scalar> potassium = {potassium_params, 19, 19, -128.71233201};

    std::vector<sae_param<scalar>> krypton_params{
        {}
    };

    sae<scalar> hydrogen = {std::vector<sae_param<scalar>>( 0 ), 1, 1, -.5};


    // call function to find all the energy states here:
    if ( params.atom() == "hydrogen" ) {
        if ( params.zeff() != 0.0 ) {
            hydrogen.N = 0;
            hydrogen.Z = params.zeff() - 1;
            hydrogen.gs_energy += std::pow( 2, params.zeff() );
        }

        hydrogen.N -= params.charge();
        if ( params.charge() > 0 ) {
            std::cerr << "hydrogen with zero electrons is just stupid..."
                      << std::endl;
            throw std::exception();
        }
        if ( params.fs() )
            numerov::find_basis_set<scalar>(
                ( memoized_finestructure_pot<scalar>( hydrogen ) ), params,
                hydrogen );
        else
            numerov::find_basis_set<scalar>(
                ( memoized_pot<scalar>( hydrogen ) ), params, hydrogen );
    }
    if ( params.atom() == "helium" ) {
        helium.N -= params.charge();
        if ( params.charge() > 0 )
            helium.gs_energy *= std::pow( 2, params.charge() );
        if ( params.fs() )
            numerov::find_basis_set<scalar>(
                ( memoized_finestructure_pot<scalar>( helium ) ), params,
                helium );
        else
            numerov::find_basis_set<scalar>( ( memoized_pot<scalar>( helium ) ),
                                             params, helium );
    }
    if ( params.atom() == "helium_tf" ) {
        numerov::find_basis_set<scalar>( &helium_tf<scalar>, params, helium );
    } else if ( params.atom() == "argon_tf" ) {
        numerov::find_basis_set<scalar>( &argon_tf<scalar>, params, argon );
    } else if ( params.atom() == "krypton" ) {
        numerov::find_basis_set<scalar>( &krypton_michelle<scalar>, params, argon );
    } else if ( params.atom() == "argon" ) {
        argon.N -= params.charge();
        if ( params.charge() > 0 )
            argon.gs_energy *= std::pow( 2, params.charge() );
        if ( params.fs() )
            numerov::find_basis_set<scalar>(
                ( memoized_finestructure_pot<scalar>( argon ) ), params,
                argon );
        else
            numerov::find_basis_set<scalar>( ( memoized_pot<scalar>( argon ) ),
                                             params, argon );
    } else if ( params.atom() == "neon" ) {
        neon.N -= params.charge();
        if ( params.charge() > 0 )
            neon.gs_energy *= std::pow( 2, params.charge() );
        if ( params.fs() )
            numerov::find_basis_set<scalar>(
                ( memoized_finestructure_pot<scalar>( neon ) ), params, neon );
        else
            numerov::find_basis_set<scalar>( ( memoized_pot<scalar>( neon ) ),
                                             params, neon );
    } else if ( params.atom() == "rubidium" ) {
        rubidium.N -= params.charge();
        if ( params.charge() > 0 )
            rubidium.gs_energy *= std::pow( 2, params.charge() );
        if ( params.fs() )
            numerov::find_basis_set<scalar>(
                ( memoized_finestructure_pot<scalar>( rubidium ) ), params,
                rubidium );
        else
            numerov::find_basis_set<scalar>(
                ( memoized_pot<scalar>( rubidium ) ), params, rubidium );
    } else if ( params.atom() == "potassium" ) {
        potassium.N -= params.charge();
        if ( params.charge() > 0 )
            potassium.gs_energy *= std::pow( 2, params.charge() );
        if ( params.fs() )
            numerov::find_basis_set<scalar>(
                ( memoized_finestructure_pot<scalar>( potassium ) ), params,
                potassium );
        else
            numerov::find_basis_set<scalar>(
                ( memoized_pot<scalar>( potassium ) ), params, potassium );
    }

    // numerov::find_basis_set<scalar>(
    // std::bind(parameterized_finestructure_pot<scalar>, std::placeholders::_1,
    // std::cref(hydrogen), std::placeholders::_2)
    //,params, hydrogen);
    // write out parameters:
    params.save_parameters();

    // delete the params...
    // delete params;
    SlepcFinalize();
    return 0;
}
