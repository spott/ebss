
// ebss:
#include <common/parameters/HamiltonianParameters.hpp>
#include <common/common.hpp>
#include <common/math.hpp>

// petsc:
#include <petsc.h>

// stl:
#include <vector>
#include <sstream>
#include <string>
#include <iostream>
#include <unordered_map>
#include <memory>
#include <functional>

// gsl:
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>


template <typename ReturnType, typename Arg>
std::function<ReturnType( Arg )>
memoize( std::function<ReturnType( Arg )> func )
{
    std::unordered_map<Arg, ReturnType> cache;
    return ( [=]( Arg t ) mutable {
        if ( cache.size() > 400 ) cache.clear();
        if ( cache.find( t ) == cache.end() )
            cache.emplace( t, func( t ) );
        return cache[t];
    } );
}

template <typename ReturnType, typename Arg>
std::function<ReturnType(Arg, bool)>
half_memoize( std::function<ReturnType( Arg )> func, size_t size )
{
    std::unordered_map<Arg, ReturnType> cache;
    return ( [=]( Arg t, bool remember ) mutable {
        if ( cache.find( t ) == cache.end() ) {
            // if (remember == false)
            // return func(t);
            if ( cache.size() > size ) {
                if ( remember ) {
                    cache.clear();
                    cache.emplace( t, func( t ) );
                } else
                    return func( t );
            }
            if ( cache.size() > size / 2 && !remember ) return func( t );
            cache.emplace( t, func( t ) );
        }
        return cache[t];
    } );
}

struct dipole_params
{
    BasisID a, b;
};

double dipole_matrix_element( double r, void* params )
{
    dipole_params* p = static_cast<dipole_params*>( params );

    // std::cout <<  p->a << " -- " << p->b << " -- " << r << std::endl;

    gsl_sf_result left;
    gsl_sf_result right;
    auto error_handler = gsl_set_error_handler_off();
    // auto left = gsl_sf_hydrogenicR (p->a.n, p->a.l, 1, r);
    // auto right = gsl_sf_hydrogenicR (p->b.n, p->b.l, 1, r);
    auto err = gsl_sf_hydrogenicR_e( p->a.n, p->a.l, 1, r, &left );
    if ( err == GSL_EUNDRFLW )
        left.val = 0;
    else if ( err != 0 ) {
        std::cout << err << " " << gsl_strerror( err ) << std::endl;
        throw std::exception();
    }

    err = 0;
    err = gsl_sf_hydrogenicR_e( p->b.n, p->b.l, 1, r, &right );
    if ( err == GSL_EUNDRFLW )
        right.val = 0;
    else if ( err != 0 ) {
        std::cout << err << " " << gsl_strerror( err ) << std::endl;
        throw std::exception();
    }

    gsl_set_error_handler( error_handler );

    return left.val * r * r * r * right.val;
}

int main( int argc, const char** argv )
{
    int ac = argc;
    char** av = new char* [argc];
    for ( size_t i = 0; i < argc; i++ ) {
        av[i] = new char[strlen( argv[i] )];
        std::copy( argv[i], argv[i] + strlen( argv[i] ), av[i] );
    }
    PetscInitialize( &ac, &av, PETSC_NULL, PETSC_NULL );

    HamiltonianParameters<PetscReal> params( argc, argv, MPI_COMM_WORLD );
    if ( params.rank() == 0 ) std::cout << params.print();

    // decltype( params.basis_parameters()->grid() ) grid;
    std::vector<double> grid;
    auto prototype = params.prototype();

    grid = params.basis_parameters()->grid();

    std::function<bool(int, int)> dipole_selection_rules;
    if ( params.fs() ) {
        dipole_selection_rules = [prototype]( int i, int j ) {
            // since j is multiplied by 2, need to have a 2 between them.
            return (
                ( std::abs( prototype[i].l - prototype[j].l ) == 1 &&
                  ( std::abs( prototype[i].j - prototype[j].j ) == 2 ||
                    prototype[i].j == prototype[j].j ) ) ||
                i == j );
        };
    } else {
        dipole_selection_rules = [prototype]( int i, int j ) {
            return ( std::abs( prototype[i].l - prototype[j].l ) == 1 ||
                     i == j );
        };
    }

    std::function<Range<typename std::vector<double>::iterator>( BasisID )>
    import_wf = [&params,
                 &grid ]( BasisID a )
                            ->Range<typename std::vector<double>::iterator>
    {
        typedef typename std::vector<double>::iterator iterator;
        static std::vector<double> l_block1;
        static int l1 = -1;
        static std::vector<double> l_block2;
        static int l2 = -1;
        static bool b = false;
        // std::cout << "l1: " << l1 << " l2: " << l2 << " l " << a.l << "
        // b " << b << std::endl;
        if ( a.l == l1 )
            return Range<iterator>{
                l_block1.begin() + (a.n - ( a.l + 1 ) ) * grid.size(),
                l_block1.begin() + (a.n - ( a.l + 1 ) ) * grid.size() +
                    grid.size()};
        else if ( a.l == l2 )
            return Range<iterator>{
                l_block2.begin() + ( a.n - ( a.l + 1 ) ) * grid.size(),
                l_block2.begin() + ( a.n - ( a.l + 1 ) ) * grid.size() +
                    grid.size()};
        else if ( a.l != l1 && a.l != l2 && b ) {
            l2 = a.l;
            l_block2 = common::import_vector_binary<double>(
                params.basis_parameters()->l_block_filename( a.l ) );
            b = false;
            return Range<iterator>{
                l_block2.begin() + ( a.n - ( a.l + 1) ) * grid.size(),
                l_block2.begin() + ( a.n - ( a.l + 1) ) * grid.size() +
                    grid.size()};
        } else if ( a.l != l1 && a.l != l2 && !b ) {
            l1 = a.l;
            l_block1 = common::import_vector_binary<double>(
                params.basis_parameters()->l_block_filename( a.l ) );
            b = true;
            return Range<iterator>{
                l_block1.begin() + ( a.n - ( a.l + 1 ) ) * grid.size(),
                l_block1.begin() + ( a.n - ( a.l + 1 ) ) * grid.size() +
                    grid.size()};
        }
    };


    std::function<PetscScalar(int, int)> findvalue;
    if ( params.fs() ) {
        findvalue = [&prototype, &grid, &import_wf ]( int i, int j )
                                                        ->PetscScalar
        {
            if ( i == j ) return 0.0;
            auto a = import_wf( prototype[i] );
            auto b = import_wf( prototype[j] );
            int s = 0;

            PetscScalar radial =
                math::integrateTrapezoidRule( a, b, grid );
            PetscScalar angular = math::CGCoefficient<PetscScalar>(
                prototype[i], prototype[j] );
            // we are only considering spin up electrons.  so m_j always ==
            // +1/2 (since m_l is always 0)
            angular *= std::sqrt(
                ( prototype[i].l +
                  ( prototype[i].j - 2 * prototype[i].l ) * .5 * .5 +
                  .5 ) *
                ( prototype[j].l +
                  ( prototype[j].j - 2 * prototype[j].l ) * .5 * .5 +
                  .5 ) /
                ( ( 2. * prototype[i].l + 1 ) *
                  ( 2. * prototype[j].l + 1 ) ) );
            if ( prototype[i].n == 2 && prototype[j].n == 2 )
                std::cout << "n = 2 transition: " << radial* angular
                          << std::endl;
            if (angular != angular || radial != radial) //check for NaNs
            {
                std::cerr << " got a NaN @ (" << i << ", " << j << "): " << prototype[i] << " <=> " << prototype[j] << std::endl;
                throw std::exception();
            }
            return radial * angular;
        };
    } else {
        findvalue = [&prototype, &grid, &import_wf ]( int i, int j )
                                                        ->PetscScalar
        {
            if ( i == j ) return 0.0;
            auto a = import_wf( prototype[i] );
            auto b = import_wf( prototype[j] );

            PetscScalar radial =
                math::integrateTrapezoidRule( a, b, grid );
            PetscScalar angular = math::CGCoefficient<PetscScalar>(
                prototype[i], prototype[j] );
            if ( prototype[i].n == 2 && prototype[j].n == 2 )
                std::cout << "n = 2 transition: " << radial* angular
                          << std::endl;
            if (angular != angular || radial != radial) //check for NaNs
            {
                std::cerr << " got a NaN @ (" << i << ", " << j << "): " << prototype[i] << " <=> " << prototype[j] << std::endl;
                throw std::exception();
            }
            return radial * angular;
        };
    }


    Mat H = common::populate_matrix<std::complex<double>>(
        params,
        dipole_selection_rules,
        findvalue,
        prototype.size(),
        prototype.size(),
        true );


    params.write_dipole_matrix( H );
    params.write_energy_eigenvalues();

    params.save_parameters();

    auto maximum_n_it = std::find_if(prototype.begin(), prototype.end(), [&params]( BasisID a ) { return a.n == params.nmax() && a.l == 0; });
    int loc = maximum_n_it - prototype.begin();
    //find the appropriate absorber size vector:
    Vec A, B;
    MatGetVecs(H, &A, &B);
    VecSetValue(A,loc,1.,INSERT_VALUES);
    VecAssemblyBegin(A);
    VecAssemblyEnd(A);

    MatMult(H, A, B);

    PetscViewer view;

    PetscViewerBinaryOpen( params.comm(),
            (params.hamiltonian_folder() + "/abs.dat").c_str(),
            FILE_MODE_WRITE,
            &view );
    VecView( B, view );

    // delete params;
    PetscFinalize();
    return 0;
}
