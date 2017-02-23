#pragma once

#include <cassert>
#include <type_traits>
#include <vector>

#include <common/common.hpp>
#include <common/parameters/BasisParameters.hpp>
#include <common/parameters/Parameters.hpp>
#include <common/types.hpp>
#include <findbasis/single_active_electron.hpp>
#include <slepc.h>


namespace finite_difference
{
template <typename scalar>
std::vector<BasisID> find_basis( std::function<scalar( scalar, BasisID )> pot,
                                 BasisParameters<scalar>& params,
                                 const sae<scalar>& atom, const BasisID& state )
{
    int order = 2;
    // dr is the number of points - the value of the grid at the beginning.
    scalar dr             = ( params.rmax() ) / ( params.points() );
    int    absorber_size  = 0;
    scalar absorber_angle = .2;

    std::vector<std::complex<scalar>> rgrid;

    rgrid.reserve( params.points() );

    size_t i = 0;
    for ( auto& r : rgrid ) {
        r = i * dr;
        if ( i > params.points() - absorber_size )
            r = ( params.points() - absorber_size ) +
                ( r - ( params.points() - absorber_size ) ) *
                    std::exp( std::complex<double>( 0, absorber_angle ) );
        i++;
    }
    std::vector<BasisID>* energies = params.basis_prototype();
    energies->resize( 0 );


    std::function<bool( int, int )> t = [order]( int i, int j ) {
        return ( std::abs( i - j ) < order || std::abs( j - i ) < order );
    };
    std::function<scalar( scalar )> potential = [pot, state]( scalar r ) {
        return pot( r, state ) + state.l * ( state.l + 1 ) / ( 2 * r * r );
    };

    std::function<scalar( int, int )> fv = [potential, dr, rgrid]( int i,
                                                                   int j ) {
        if ( i == j )  // diagonal
            return potential( rgrid->at( i ) ) + 1 / ( dr * dr );
        else if ( i - 1 == j || i + 1 == j )  // off diagonal
            return -1 / ( 2 * dr * dr );
    };

    Mat H =
        common::populate_matrix( static_cast<const Parameters>( params ), t, fv,
                                 params.points(), params.points(), true );

    Vec xr, xi;
    MatGetVecs( H, PETSC_NULL, &xr );
    MatGetVecs( H, PETSC_NULL, &xi );

    EPS eps;
    EPSCreate( params.comm(), &eps );

    EPSSetOperators( eps, H, PETSC_NULL );
    EPSSetProblemType( eps, EPS_HEP );
    EPSSetDimensions( eps, params.nmax() - state.l, PETSC_DECIDE,
                      PETSC_DECIDE );
    EPSSetWhichEigenpairs( eps, EPS_SMALLEST_REAL );
    EPSSetFromOptions( eps );
    if ( params.rank() == 0 ) std::cout << "starting solve" << std::endl;
    EPSSolve( eps );

    PetscInt    its, nconv, maxit, nev;
    PetscReal   error, tol, re, im;
    PetscScalar kr, ki;
    EPSGetIterationNumber( eps, &its );
    PetscPrintf( params.comm(), " Number of iterations of the method: %D\n",
                 its );
    EPSGetDimensions( eps, &nev, PETSC_NULL, PETSC_NULL );
    PetscPrintf( params.comm(), " Number of requested eigenvalues: %D\n", nev );
    EPSGetTolerances( eps, &tol, &maxit );
    PetscPrintf( params.comm(), " Stopping condition: tol=%.4G, maxit=%D\n",
                 tol, maxit );

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Display solution and clean up
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       */
    /*
       Get number of converged approximate eigenpairs
       */
    EPSGetConverged( eps, &nconv );
    PetscPrintf( params.comm(), " Number of converged eigenpairs: %D\n\n",
                 nconv );

    BasisID tmp = state;

    if ( nconv > 0 ) {
        /*
           Display eigenvalues and relative errors
           */
        PetscPrintf( params.comm(),
                     "           k          ||Ax-kx||/||kx||\n"
                     "   ----------------- ------------------\n" );

        for ( int i = 0; i < nconv; i++ ) {
            /*
               Get converged eigenpairs: i-th eigenvalue is stored in kr
               (real part) and
               ki (imaginary part)
               */
            EPSGetEigenpair( eps, i, &kr, &ki, xr, xi );
            /*
               Compute the relative error associated to each eigenpair
               */
            EPSComputeRelativeError( eps, i, &error );

            re = PetscRealPart( kr );
            im = PetscImaginaryPart( kr );


            if ( im != 0.0 ) {
                PetscPrintf( params.comm(), " %9F%+9F j %12G\n", re, im,
                             error );
            } else {
                PetscPrintf( params.comm(), "   %12F       %12G\n", re, error );
            }
            tmp.e = kr;
            tmp.n = i + 1 + tmp.l;
            energies->push_back( tmp );
            auto vout1 = common::Vec_to_vector( xr );
            if ( params.rank() == 0 )
                common::export_vector_binary<PetscReal>(
                    params.basis_function_filename( tmp ),
                    common::vector_type_change<PetscScalar, PetscReal>(
                        vout1 ) );
        }
        PetscPrintf( params.comm(), "\n" );
    }
    MatDestroy( &H );
    EPSDestroy( &eps );
    VecDestroy( &xr );
    VecDestroy( &xi );
}


// template <typename scalar, typename write_type>
// void find_basis_set( std::function<scalar( scalar, BasisID )> potential,
//                      BasisParameters<scalar, write_type>& params,
//                      sae<scalar> atom )
// {
//     // auto rgrid =
//     std::vector<BasisID> energies;

//     for ( size_t l = 0; l <= params.lmax(); ++l ) {
//         for ( int j = ( ( l > 0 ) ? 2 * l - 1 : 1 );
//               j <= ( ( l > 0 ) ? 2 * l + 1 : 1 );
//               j += 2 ) {
//             tmp.n = 0;
//             tmp.l = l;
//             tmp.j = j;
//             tmp.m = 0;
//             tmp.e = 0;
//             e_tmp = find_basis( potential, params, atom, tmp );
//         }
//     }
// }
}
