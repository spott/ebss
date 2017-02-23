
#include <petsc.h>

// stl:
#include <cassert>
#include <complex>
#include <utility>
#include <vector>

// mine:
#include <common/common.hpp>
#include <common/math.hpp>
#include <common/parameters/HamiltonianParameters.hpp>
#include <common/parameters/NonlinearParameters.hpp>
#include <common/parameters/StateParameters.hpp>


Vec psi( int order, std::vector<double>::const_iterator frequencies_begin,
         std::vector<double>::const_iterator frequencies_end, PetscScalar wg,
         Vec& H0, Mat& D, Vec& psi0, Vec& mask,
         std::vector<BasisID>& prototype );
Vec psi_conjugate( int                                 order,
                   std::vector<double>::const_iterator frequencies_begin,
                   std::vector<double>::const_iterator frequencies_end,
                   PetscScalar wg, Vec& H0, Mat& D, Vec& psi0, Vec& mask,
                   std::vector<BasisID>& prototype );

struct VectorElement {
    const std::tuple<std::complex<double>, int>& el;
    const std::vector<BasisID>& proto;

    VectorElement( const std::tuple<std::complex<double>, int>& element,
                   const std::vector<BasisID>& prototype )
        : el( element ), proto( prototype ){};
};

std::ostream& operator<<( std::ostream& out, const VectorElement a )
{
    out << std::get<0>( a.el ) << "\t" << a.proto[std::get<1>( a.el )].n << "\t"
        << a.proto[std::get<1>( a.el )].l;
    return out;
}

namespace std
{
bool operator<( const std::tuple<std::complex<double>, int>& a,
                const std::tuple<std::complex<double>, int>& b )
{
    return std::abs( std::get<0>( a ).real() ) <
           std::abs( std::get<0>( b ).real() );
}
}

int main( int argc, const char** argv )
{
    int    ac = argc;
    char** av = new char*[argc + 1];
    for ( int i = 0; i < argc; i++ ) {
        av[i] = new char[strlen( argv[i] ) + 1];
        std::copy( argv[i], argv[i] + strlen( argv[i] ) + 1, av[i] );
    }
    av[argc] = NULL;
    PetscInitialize( &ac, &av, PETSC_NULL, PETSC_NULL );

    // the parameters (where to find the hamiltonian)
    NonlinearParameters           nparams( argc, argv, PETSC_COMM_WORLD );
    HamiltonianParameters<double> params( argc, argv, PETSC_COMM_WORLD );
    StateParameters               sparams( argc, argv, PETSC_COMM_WORLD );

    auto prototype = params.prototype();

    MPI_Barrier( PETSC_COMM_WORLD );

    std::cout << std::scientific;

    if ( params.rank() == 0 ) {
        std::
        cout << "#Git commit: " << GIT_COMMIT << std::endl;
        std::cout << "reading in matrix:" << std::endl;
    }

    // read in the matrices
    Mat D = params.read_dipole_matrix();
    MatAssemblyBegin( D, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd( D, MAT_FINAL_ASSEMBLY );
    if ( params.rank() == 0 ) std::cout << "reading in vector:" << std::endl;
    Vec H0 = params.read_energy_eigenvalues();
    VecAssemblyBegin( H0 );
    VecAssemblyEnd( H0 );
    int size;
    VecGetSize( H0, &size );
    // std::cout << size << std::endl;

    assert( size == int( prototype.size() ) );
    std::function<bool( int, int )> dipole_selection_rules = [prototype](
        int i, int j ) {
        return ( std::abs( prototype[i].l - prototype[j].l ) == 1 || i == j );
    };

    // Mat D = common::populate_matrix< std::complex<double> >(params,
    // dipole_selection_rules,
    //[](int i, int j) {return 1.; },
    // prototype.size(),
    // prototype.size(), true);

    // Vec maximums;
    // VecDuplicate(H0, &maximums);
    // int* locations;
    // MatGetRowMax(D,maximums,PETSC_NULL);
    // VecView(maximums, PETSC_VIEWER_STDOUT_WORLD);

    // create a mask for bound states:


    // Vec mask = common::map_function(H0, [](PetscScalar in) { return 1; });

    auto imgs = nparams.imgs();
    Vec  psi0;
    Vec  psi1;

    std::vector<std::tuple<PetscScalar, int>> maxes;
    // the final "psi"
    if ( nparams.wf_filename().empty() ) {
        VecDuplicate( H0, &psi1 );
        VecSetValue( psi1, 0, 1., INSERT_VALUES );
        VecAssemblyBegin( psi1 );
        VecAssemblyEnd( psi1 );
        maxes.push_back( std::make_tuple( std::complex<double>{1, 0}, 0 ) );
    } else {
        psi1 = common::petsc_binary_read<Vec>( nparams.wf_filename(),
                                               params.comm() );
        // Sort vector:
        maxes = math::VecStupidSort( psi1, []( PetscScalar a, PetscScalar b ) {
            return std::abs( a ) > std::abs( b );
        } );
    }


    if ( params.rank() == 0 )
        for ( auto b : maxes ) {
            std::cout << std::get<0>( b ) << ", " << std::get<1>( b )
                      << std::endl;
        }

    // the mask
    Vec mask;
    VecDuplicate( psi1, &mask );

    // the initial psi
    VecDuplicate( H0, &psi0 );
    VecSetValue( psi0, 0, 1., INSERT_VALUES );
    VecAssemblyBegin( psi0 );
    VecAssemblyEnd( psi0 );

    // find the "ground state":
    VecCopy( psi0, mask );
    PetscScalar wg;
    VecPointwiseMult( mask, psi0, H0 );
    VecDot( mask, psi0, &wg );

    // the mask should be a copy of the starting psi
    VecCopy( psi0, mask );
    // get the frequencies list from nonlinear params:
    auto freqs = nparams.freqs();
    std::sort( freqs.begin(), freqs.end() );

    try {
        common::export_vector_binary( std::string( "./frequencies.dat" ),
                                      freqs );
    } catch ( ... ) {
        std::cout << "error in export_vector_binary" << std::endl;
    }


    std::vector<std::vector<std::complex<double>>> chi1_data(
        nparams.chi1s().size() );
    for ( auto a : chi1_data ) a.reserve( freqs.size() * imgs.size() );

    std::vector<std::vector<std::complex<double>>> chi3_data(
        nparams.chi3s().size() );
    for ( auto a : chi3_data ) a.reserve( freqs.size() * imgs.size() );

    std::vector<std::vector<std::complex<double>>> chi5_data(
        nparams.chi5s().size() );
    for ( auto a : chi5_data ) a.reserve( freqs.size() * imgs.size() );

    std::vector<std::vector<std::complex<double>>> chi7_data(
        nparams.chi7s().size() );
    for ( auto a : chi7_data ) a.reserve( freqs.size() * imgs.size() );

    std::vector<std::vector<std::complex<double>>> chi9_data(
        nparams.chi9s().size() );
    for ( auto a : chi9_data ) a.reserve( freqs.size() * imgs.size() );

    std::vector<std::vector<std::complex<double>>> chi11_data(
        nparams.chi11s().size() );
    for ( auto a : chi11_data ) a.reserve( freqs.size() * imgs.size() );

    for ( auto j : imgs ) {
        VecShift( H0, std::complex<double>( 0, -j ) );
        PetscScalar wg;
        VecDot( psi0, H0, &wg );
        if ( params.rank() == 0 ) std::cout << j << " wg: " << wg << std::endl;

        for ( size_t f = 0; f < freqs.size(); ++f ) {
            if ( params.rank() == 0 )
                std::cout << f << "(" << freqs[f] << ")" << std::endl;
            PetscScalar t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
            Vec         c;
            VecDuplicate( H0, &c );
            // we need to loop through all permutations of all chi types:

            // first we need to generate the correct vectors:
            // Chi1
            for ( auto i = nparams.chi1s().begin(); i != nparams.chi1s().end();
                  ++i ) {
                PetscScalar total = 0;
                for ( auto max1 : maxes ) {
                    for ( auto max2 : maxes ) {
                        VecSet( psi0, 0 );
                        VecSet( psi1, 0 );

                        VecSetValue( psi0, std::get<1>( max2 ), 1.,
                                     INSERT_VALUES );
                        VecSetValue( psi1, std::get<1>( max1 ), 1.,
                                     INSERT_VALUES );
                        VecAssemblyBegin( psi0 );
                        VecAssemblyBegin( psi1 );
                        VecAssemblyEnd( psi0 );
                        VecAssemblyEnd( psi1 );
                        PetscScalar wg0;
                        VecPointwiseMult( mask, psi0, H0 );
                        VecDot( mask, psi0, &wg0 );

                        PetscScalar wg1;
                        VecPointwiseMult( mask, psi1, H0 );
                        VecDot( mask, psi1, &wg1 );

                        if ( params.rank() == 0 )
                            std::cout
                                << "========================================="
                                << std::endl
                                << *i << ": (" << std::get<1>( max1 ) << ","
                                << std::get<1>( max2 ) << ") wg: " << wg0
                                << ", " << wg1 << std::endl;
                        std::vector<double> freq{( *i ) * freqs[f]};
                        Vec p0 = psi( 0, freq.cbegin(), freq.cend(), wg0, H0, D,
                                      psi0, psi0, prototype );
                        Vec p1 = psi( 1, freq.cbegin(), freq.cend(), wg0, H0, D,
                                      psi0, psi0, prototype );
                        Vec p1c =
                            psi_conjugate( 1, freq.cbegin(), freq.cend(), wg1,
                                           H0, D, psi1, psi1, prototype );
                        Vec p0c =
                            psi_conjugate( 0, freq.cbegin(), freq.cend(), wg1,
                                           H0, D, psi1, psi1, prototype );

                        // chi1 = <\psi^(0) | D | \psi^(1)>
                        MatMult( D, p1, c );
                        VecDot( p0c, c, &t1 );
                        // chi1 = <\psi^(1) | D | \psi^(0)>
                        MatMult( D, p0, c );
                        VecDot( p1c, c, &t2 );
                        total += ( t1 + t2 ) *
                                 std::conj( std::get<0>( max1 ) ) *
                                 std::get<0>( max2 );
                        if ( params.rank() == 0 &&
                             std::abs( t1 + t2 ) >= 1e-16 )
                            std::cout
                                << "========================================="
                                << std::endl
                                << *i << ": (" << std::get<1>( max1 ) << ","
                                << std::get<1>( max2 ) << ") wg: " << wg0
                                << ", " << wg1 << std::endl;
                        if ( params.rank() == 0 &&
                             std::abs( t1 + t2 ) >= 1e-16 )
                            std::cout << "terms ("
                                      << prototype[std::get<1>( max1 )] << "->"
                                      << prototype[std::get<1>( max2 )]
                                      << "): " << t1 << ", " << t2 << std::endl;
                        if ( params.rank() == 0 &&
                             std::abs( t1 + t2 ) >= 1e-16 )
                            std::cout << "final: ("
                                      << prototype[std::get<1>( max1 )] << "->"
                                      << prototype[std::get<1>( max2 )] << "): "
                                      << ( t1 + t2 ) *
                                             std::conj( std::get<0>( max1 ) ) *
                                             std::get<0>( max2 )
                                      << std::endl;
                        VecDestroy( &p1 );
                        VecDestroy( &p0 );
                        VecDestroy( &p1c );
                        VecDestroy( &p0c );
                    }
                }
                if ( params.rank() == 0 )
                    std::cout << "final: " << ( total ) << std::endl;
                chi1_data[i - nparams.chi1s().begin()].push_back( ( total ) );
            }

            // Chi3
            for ( auto i = nparams.chi3s().begin(); i != nparams.chi3s().end();
                  ++i ) {
                PetscScalar total = 0;
                for ( auto max1 : maxes ) {
                    for ( auto max2 : maxes ) {
                        if ( std::abs( prototype[std::get<1>( max1 )].l -
                                       prototype[std::get<1>( max2 )].l ) %
                                 2 !=
                             0 )
                            break;
                        VecSet( psi0, 0 );
                        VecSet( psi1, 0 );

                        VecSetValue( psi0, std::get<1>( max2 ), 1.,
                                     INSERT_VALUES );
                        VecSetValue( psi1, std::get<1>( max1 ), 1.,
                                     INSERT_VALUES );
                        VecAssemblyBegin( psi0 );
                        VecAssemblyBegin( psi1 );
                        VecAssemblyEnd( psi0 );
                        VecAssemblyEnd( psi1 );
                        PetscScalar wg0;
                        VecPointwiseMult( mask, psi0, H0 );
                        VecDot( mask, psi0, &wg0 );

                        PetscScalar wg1;
                        VecPointwiseMult( mask, psi1, H0 );
                        VecDot( mask, psi1, &wg1 );
                        if ( params.rank() == 0 )
                            std::cout
                                << "========================================="
                                << std::endl
                                << *i << ": (" << prototype[std::get<1>( max1 )]
                                << " [" << prototype[std::get<0>( max1 )]
                                << "]," << std::get<1>( max2 ) << " ["
                                << std::get<0>( max2 ) << "]) wg: " << wg0
                                << ", " << wg1 << std::endl;
                        std::sort( ( *i ).begin(), ( *i ).end() );
                        size_t multiplicity = 1;
                        std::array<int, 3> ts{{0, 0, 0}};
                        for ( auto m : ( *i ) ) {
                            if ( m == -1 ) ts[0]++;
                            if ( m == 1 ) ts[2]++;
                            if ( m == 0 ) ts[1]++;
                        }
                        for ( auto m : ts )
                            multiplicity *= math::factorial( m );

                        PetscScalar result = 0;
                        do {
                            std::vector<double> freq{( *i )[0] * freqs[f],
                                                     ( *i )[1] * freqs[f],
                                                     ( *i )[2] * freqs[f]};
                            Vec p3 = psi( 3, freq.cbegin(), freq.cend(), wg0,
                                          H0, D, psi0, psi0, prototype );
                            Vec p2 = psi( 2, freq.cbegin(), freq.cend(), wg0,
                                          H0, D, psi0, psi0, prototype );
                            Vec p0 = psi( 0, freq.cbegin(), freq.cend(), wg0,
                                          H0, D, psi0, psi0, prototype );
                            Vec p1 = psi( 1, freq.cbegin(), freq.cend(), wg0,
                                          H0, D, psi0, psi0, prototype );
                            Vec p3c = psi_conjugate( 3, freq.cbegin(),
                                                     freq.cend(), wg1, H0, D,
                                                     psi1, psi1, prototype );
                            Vec p2c = psi_conjugate( 2, freq.cbegin(),
                                                     freq.cend(), wg1, H0, D,
                                                     psi1, psi1, prototype );
                            Vec p1c = psi_conjugate( 1, freq.cbegin(),
                                                     freq.cend(), wg1, H0, D,
                                                     psi1, psi1, prototype );
                            Vec p0c = psi_conjugate( 0, freq.cbegin(),
                                                     freq.cend(), wg1, H0, D,
                                                     psi1, psi1, prototype );

                            // chi3 = <\psi1^(0) | D | \psi0^(3)>
                            MatMult( D, p3, c );
                            VecDot( p0c, c, &t1 );
                            //      + <\psi1^(1) | D | \psi0^(2)>
                            MatMult( D, p2, c );
                            VecDot( p1c, c, &t2 );
                            //      + <\psi1^(2) | D | \psi0^(1)>
                            MatMult( D, p1, c );
                            VecDot( p2c, c, &t3 );
                            //      + <\psi1^(3) | D | \psi0^(0)>
                            MatMult( D, p0, c );
                            VecDot( p3c, c, &t4 );
                            result += ( t1 + t2 + t3 + t4 ) *
                                      std::conj( std::get<0>( max1 ) ) *
                                      std::get<0>( max2 );
                            if ( params.rank() == 0 &&
                                 std::abs( t1 + t2 + t3 + t4 ) >= 1e-16 )
                                std::cout
                                    << "terms: ("
                                    << prototype[std::get<1>( max1 )] << ","
                                    << prototype[std::get<1>( max2 )]
                                    << "): " << t1 << ", " << t2 << ", " << t3
                                    << ", " << t4 << " || "
                                    << std::conj( std::get<0>( max1 ) ) *
                                           std::get<0>( max2 )
                                    << std::endl;
                            VecDestroy( &p3 );
                            VecDestroy( &p2 );
                            VecDestroy( &p1 );
                            VecDestroy( &p0 );
                            VecDestroy( &p3c );
                            VecDestroy( &p2c );
                            VecDestroy( &p1c );
                            VecDestroy( &p0c );
                        } while ( std::next_permutation( ( *i ).begin(),
                                                         ( *i ).end() ) );
                        if ( params.rank() == 0 &&
                             std::abs( t1 + t2 + t3 + t4 ) >= 1e-16 )
                            std::cout << "final: ("
                                      << prototype[std::get<1>( max1 )] << ","
                                      << prototype[std::get<1>( max2 )]
                                      << "): " << result << " multiplicity "
                                      << multiplicity << std::endl;
                        total += result * static_cast<double>( multiplicity ) /
                                 static_cast<double>( math::factorial( 3 ) );
                    }
                }
                if ( params.rank() == 0 )
                    std::cout << "final: " << total << std::endl;
                chi3_data[i - nparams.chi3s().begin()].push_back( total );
            }
            // chi5
            for ( auto i = nparams.chi5s().begin(); i != nparams.chi5s().end();
                  ++i ) {
                if ( params.rank() == 0 )
                    std::cout << "========================================="
                              << std::endl
                              << *i << std::endl;
                std::sort( ( *i ).begin(), ( *i ).end() );
                size_t multiplicity = 1;
                std::array<int, 3> ts{{0, 0, 0}};
                for ( auto m : ( *i ) ) {
                    if ( m == -1 ) ts[0]++;
                    if ( m == 0 ) ts[1]++;
                    if ( m == 1 ) ts[2]++;
                }
                for ( auto  m : ts ) multiplicity *= math::factorial( m );
                PetscScalar result;
                do {
                    std::vector<double> freq{
                        ( *i )[0] * freqs[f], ( *i )[1] * freqs[f],
                        ( *i )[2] * freqs[f], ( *i )[3] * freqs[f],
                        ( *i )[4] * freqs[f]};
                    // if (params.rank() == 0) std::cout << "p5: " << std::endl;
                    Vec p5 = psi( 5, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    // if (params.rank() == 0) std::cout << "p4: " << std::endl;
                    Vec p4 = psi( 4, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    // if (params.rank() == 0) std::cout << "p3: " << std::endl;
                    Vec p3 = psi( 3, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    // if (params.rank() == 0) std::cout << "p2: " << std::endl;
                    Vec p2 = psi( 2, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    // if (params.rank() == 0) std::cout << "p1: " << std::endl;
                    Vec p1 = psi( 1, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    // if (params.rank() == 0) std::cout << "p0: " << std::endl;
                    Vec p0 = psi( 0, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    // if (params.rank() == 0) std::cout << "p5c: " <<
                    // std::endl;
                    Vec p5c = psi_conjugate( 5, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi1, mask, prototype );
                    // if (params.rank() == 0) std::cout << "p4c: " <<
                    // std::endl;
                    Vec p4c = psi_conjugate( 4, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi1, mask, prototype );
                    // if (params.rank() == 0) std::cout << "p3c: " <<
                    // std::endl;
                    Vec p3c = psi_conjugate( 3, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi1, mask, prototype );
                    // if (params.rank() == 0) std::cout << "p2c: " <<
                    // std::endl;
                    Vec p2c = psi_conjugate( 2, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi1, mask, prototype );
                    // if (params.rank() == 0) std::cout << "p1c: " <<
                    // std::endl;
                    Vec p1c = psi_conjugate( 1, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi1, mask, prototype );
                    // if (params.rank() == 0) std::cout << "p0c: " <<
                    // std::endl;
                    Vec p0c = psi_conjugate( 0, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi1, mask, prototype );

                    // chi5 = <\psi^(0) | D | \psi^(5)>
                    MatMult( D, p5, c );
                    VecDot( p0c, c, &t1 );
                    //      + <\psi^(1) | D | \psi^(4)>
                    MatMult( D, p4, c );
                    VecDot( p1c, c, &t2 );
                    //      + <\psi^(2) | D | \psi^(3)>
                    MatMult( D, p3, c );
                    VecDot( p2c, c, &t3 );
                    //      + <\psi^(2) | D | \psi^(3)>
                    MatMult( D, p2, c );
                    VecDot( p3c, c, &t4 );
                    //      + <\psi^(4) | D | \psi^(1)>
                    MatMult( D, p1, c );
                    VecDot( p4c, c, &t5 );
                    //      + <\psi^(5) | D | \psi^(0)>
                    MatMult( D, p0, c );
                    VecDot( p5c, c, &t6 );
                    if ( params.rank() == 0 )
                        std::cout << "terms: " << t1 << ", " << t2 << ", " << t3
                                  << ", " << t4 << ", " << t5 << ", " << t6
                                  << std::endl;
                    result += ( t1 + t2 + t3 + t4 );
                    VecDestroy( &p5 );
                    VecDestroy( &p4 );
                    VecDestroy( &p3 );
                    VecDestroy( &p2 );
                    VecDestroy( &p1 );
                    VecDestroy( &p0 );
                    VecDestroy( &p5c );
                    VecDestroy( &p4c );
                    VecDestroy( &p3c );
                    VecDestroy( &p2c );
                    VecDestroy( &p1c );
                    VecDestroy( &p0c );
                } while (
                    std::next_permutation( ( *i ).begin(), ( *i ).end() ) );

                if ( params.rank() == 0 )
                    std::cout << "final: " << result << std::endl;
                chi5_data[i - nparams.chi5s().begin()].push_back(
                    result * static_cast<double>( multiplicity ) /
                    static_cast<double>( math::factorial( 5 ) ) );
            }
            // Chi7
            for ( auto i = nparams.chi7s().begin(); i != nparams.chi7s().end();
                  ++i ) {
                if ( params.rank() == 0 )
                    std::cout << "========================================="
                              << std::endl
                              << *i << std::endl;
                std::sort( ( *i ).begin(), ( *i ).end() );
                size_t multiplicity = 1;
                std::array<int, 3> ts{{0, 0, 0}};
                for ( auto m : ( *i ) ) {
                    if ( m == -1 ) ts[0]++;
                    if ( m == 0 ) ts[1]++;
                    if ( m == 1 ) ts[2]++;
                }
                for ( auto  m : ts ) multiplicity *= math::factorial( m );
                PetscScalar result;
                do {
                    std::vector<double> freq{
                        ( *i )[0] * freqs[f], ( *i )[1] * freqs[f],
                        ( *i )[2] * freqs[f], ( *i )[3] * freqs[f],
                        ( *i )[4] * freqs[f], ( *i )[5] * freqs[f],
                        ( *i )[6] * freqs[f]};
                    Vec p7 = psi( 7, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p6 = psi( 6, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p5 = psi( 5, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p4 = psi( 4, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p3 = psi( 3, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p2 = psi( 2, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p1 = psi( 1, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p0 = psi( 0, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p7c = psi_conjugate( 7, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p6c = psi_conjugate( 6, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p5c = psi_conjugate( 5, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p4c = psi_conjugate( 4, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p3c = psi_conjugate( 3, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p2c = psi_conjugate( 2, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p1c = psi_conjugate( 1, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p0c = psi_conjugate( 0, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );

                    // chi7 = <\psi^(0) | D | \psi^(7)>
                    MatMult( D, p7, c );
                    VecDot( p0c, c, &t1 );
                    //      + <\psi^(1) | D | \psi^(6)>
                    MatMult( D, p6, c );
                    VecDot( p1c, c, &t2 );
                    //      + <\psi^(2) | D | \psi^(5)>
                    MatMult( D, p5, c );
                    VecDot( p2c, c, &t3 );
                    //      + <\psi^(3) | D | \psi^(4)>
                    MatMult( D, p4, c );
                    VecDot( p3c, c, &t4 );
                    //      + <\psi^(4) | D | \psi^(3)>
                    MatMult( D, p3, c );
                    VecDot( p4c, c, &t5 );
                    //      + <\psi^(5) | D | \psi^(2)>
                    MatMult( D, p2, c );
                    VecDot( p5c, c, &t6 );
                    //      + <\psi^(6) | D | \psi^(1)>
                    MatMult( D, p1, c );
                    VecDot( p6c, c, &t7 );
                    //      + <\psi^(7) | D | \psi^(0)>
                    MatMult( D, p0, c );
                    VecDot( p7c, c, &t8 );
                    if ( params.rank() == 0 )
                        std::cout << "terms: " << t1 << ", " << t2 << ", " << t3
                                  << ", " << t4 << ", " << t5 << ", " << t6
                                  << ", " << t7 << ", " << t8 << std::endl;
                    result += ( t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8 );
                    VecDestroy( &p7 );
                    VecDestroy( &p6 );
                    VecDestroy( &p5 );
                    VecDestroy( &p4 );
                    VecDestroy( &p3 );
                    VecDestroy( &p2 );
                    VecDestroy( &p1 );
                    VecDestroy( &p0 );
                    VecDestroy( &p7c );
                    VecDestroy( &p6c );
                    VecDestroy( &p5c );
                    VecDestroy( &p4c );
                    VecDestroy( &p3c );
                    VecDestroy( &p2c );
                    VecDestroy( &p1c );
                    VecDestroy( &p0c );
                } while (
                    std::next_permutation( ( *i ).begin(), ( *i ).end() ) );

                if ( params.rank() == 0 )
                    std::cout << "final: " << result << std::endl;
                chi7_data[i - nparams.chi7s().begin()].push_back(
                    result * static_cast<double>( multiplicity ) /
                    static_cast<double>( math::factorial( 7 ) ) );
            }

            // Chi9
            for ( auto i = nparams.chi9s().begin(); i != nparams.chi9s().end();
                  ++i ) {
                if ( params.rank() == 0 )
                    std::cout << "========================================="
                              << std::endl
                              << *i << std::endl;
                std::sort( ( *i ).begin(), ( *i ).end() );
                size_t multiplicity = 1;
                std::array<int, 3> ts{{0, 0, 0}};
                for ( auto m : ( *i ) ) {
                    if ( m == -1 ) ts[0]++;
                    if ( m == 0 ) ts[1]++;
                    if ( m == 1 ) ts[2]++;
                }
                for ( auto  m : ts ) multiplicity *= math::factorial( m );
                PetscScalar result;
                do {
                    std::vector<double> freq{
                        ( *i )[0] * freqs[f], ( *i )[1] * freqs[f],
                        ( *i )[2] * freqs[f], ( *i )[3] * freqs[f],
                        ( *i )[4] * freqs[f], ( *i )[5] * freqs[f],
                        ( *i )[6] * freqs[f], ( *i )[7] * freqs[f],
                        ( *i )[8] * freqs[f]};
                    Vec p9 = psi( 9, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p8 = psi( 8, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p7 = psi( 7, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p6 = psi( 6, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p5 = psi( 5, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p4 = psi( 4, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p3 = psi( 3, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p2 = psi( 2, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p1 = psi( 1, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p0 = psi( 0, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p9c = psi_conjugate( 9, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p8c = psi_conjugate( 8, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p7c = psi_conjugate( 7, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p6c = psi_conjugate( 6, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p5c = psi_conjugate( 5, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p4c = psi_conjugate( 4, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p3c = psi_conjugate( 3, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p2c = psi_conjugate( 2, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p1c = psi_conjugate( 1, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p0c = psi_conjugate( 0, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );

                    // chi9 = <\psi^(0) | D | \psi^(9)>
                    MatMult( D, p9, c );
                    VecDot( p0c, c, &t1 );
                    //      + <\psi^(1) | D | \psi^(8)>
                    MatMult( D, p8, c );
                    VecDot( p1c, c, &t2 );
                    //      + <\psi^(2) | D | \psi^(7)>
                    MatMult( D, p7, c );
                    VecDot( p2c, c, &t3 );
                    //      + <\psi^(3) | D | \psi^(6)>
                    MatMult( D, p6, c );
                    VecDot( p3c, c, &t4 );
                    //      + <\psi^(4) | D | \psi^(5)>
                    MatMult( D, p5, c );
                    VecDot( p4c, c, &t5 );
                    //      + <\psi^(5) | D | \psi^(4)>
                    MatMult( D, p4, c );
                    VecDot( p5c, c, &t6 );
                    //      + <\psi^(6) | D | \psi^(3)>
                    MatMult( D, p3, c );
                    VecDot( p6c, c, &t7 );
                    //      + <\psi^(7) | D | \psi^(2)>
                    MatMult( D, p2, c );
                    VecDot( p7c, c, &t8 );
                    //      + <\psi^(8) | D | \psi^(1)>
                    MatMult( D, p1, c );
                    VecDot( p8c, c, &t9 );
                    //      + <\psi^(9) | D | \psi^(0)>
                    MatMult( D, p0, c );
                    VecDot( p9c, c, &t10 );
                    if ( params.rank() == 0 )
                        std::cout << "terms: " << t1 << ", " << t2 << ", " << t3
                                  << ", " << t4 << ", " << t5 << ", " << t6
                                  << ", " << t7 << ", " << t8 << ", " << t9
                                  << ", " << t10 << std::endl;
                    result +=
                        ( t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9 + t10 );
                    VecDestroy( &p9 );
                    VecDestroy( &p8 );
                    VecDestroy( &p7 );
                    VecDestroy( &p6 );
                    VecDestroy( &p5 );
                    VecDestroy( &p4 );
                    VecDestroy( &p3 );
                    VecDestroy( &p2 );
                    VecDestroy( &p1 );
                    VecDestroy( &p0 );
                    VecDestroy( &p9c );
                    VecDestroy( &p8c );
                    VecDestroy( &p7c );
                    VecDestroy( &p6c );
                    VecDestroy( &p5c );
                    VecDestroy( &p4c );
                    VecDestroy( &p3c );
                    VecDestroy( &p2c );
                    VecDestroy( &p1c );
                    VecDestroy( &p0c );
                } while (
                    std::next_permutation( ( *i ).begin(), ( *i ).end() ) );

                if ( params.rank() == 0 )
                    std::cout << "final: " << result << std::endl;
                chi9_data[i - nparams.chi9s().begin()].push_back(
                    result * static_cast<double>( multiplicity ) /
                    static_cast<double>( math::factorial( 9 ) ) );
            }

            // Chi11
            for ( auto i = nparams.chi11s().begin();
                  i != nparams.chi11s().end(); ++i ) {
                if ( params.rank() == 0 )
                    std::cout << "========================================="
                              << std::endl
                              << *i << std::endl;
                std::sort( ( *i ).begin(), ( *i ).end() );
                size_t multiplicity = 1;
                std::array<int, 3> ts{{0, 0, 0}};
                for ( auto m : ( *i ) ) {
                    if ( m == -1 ) ts[0]++;
                    if ( m == 0 ) ts[1]++;
                    if ( m == 1 ) ts[2]++;
                }
                for ( auto  m : ts ) multiplicity *= math::factorial( m );
                PetscScalar result;
                do {
                    std::vector<double> freq{
                        ( *i )[0] * freqs[f], ( *i )[1] * freqs[f],
                        ( *i )[2] * freqs[f], ( *i )[3] * freqs[f],
                        ( *i )[4] * freqs[f], ( *i )[5] * freqs[f],
                        ( *i )[6] * freqs[f], ( *i )[7] * freqs[f],
                        ( *i )[8] * freqs[f], ( *i )[9] * freqs[f],
                        ( *i )[10] * freqs[f]};

                    Vec p11 = psi( 11, freq.cbegin(), freq.cend(), wg, H0, D,
                                   psi0, mask, prototype );
                    Vec p10 = psi( 10, freq.cbegin(), freq.cend(), wg, H0, D,
                                   psi0, mask, prototype );
                    Vec p9 = psi( 9, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p8 = psi( 8, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p7 = psi( 7, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p6 = psi( 6, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p5 = psi( 5, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p4 = psi( 4, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p3 = psi( 3, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p2 = psi( 2, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p1 = psi( 1, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p0 = psi( 0, freq.cbegin(), freq.cend(), wg, H0, D,
                                  psi0, mask, prototype );
                    Vec p11c =
                        psi_conjugate( 11, freq.cbegin(), freq.cend(), wg, H0,
                                       D, psi0, mask, prototype );
                    Vec p10c =
                        psi_conjugate( 10, freq.cbegin(), freq.cend(), wg, H0,
                                       D, psi0, mask, prototype );
                    Vec p9c = psi_conjugate( 9, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p8c = psi_conjugate( 8, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p7c = psi_conjugate( 7, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p6c = psi_conjugate( 6, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p5c = psi_conjugate( 5, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p4c = psi_conjugate( 4, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p3c = psi_conjugate( 3, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p2c = psi_conjugate( 2, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p1c = psi_conjugate( 1, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );
                    Vec p0c = psi_conjugate( 0, freq.cbegin(), freq.cend(), wg,
                                             H0, D, psi0, mask, prototype );

                    // chi9 = <\psi^(0) | D | \psi^(11)>
                    MatMult( D, p11, c );
                    VecDot( p0c, c, &t1 );
                    //      + <\psi^(1) | D | \psi^(10)>
                    MatMult( D, p10, c );
                    VecDot( p1c, c, &t2 );
                    //      + <\psi^(2) | D | \psi^(9)>
                    MatMult( D, p9, c );
                    VecDot( p2c, c, &t3 );
                    //      + <\psi^(3) | D | \psi^(8)>
                    MatMult( D, p8, c );
                    VecDot( p3c, c, &t4 );
                    //      + <\psi^(4) | D | \psi^(7)>
                    MatMult( D, p7, c );
                    VecDot( p4c, c, &t5 );
                    //      + <\psi^(5) | D | \psi^(6)>
                    MatMult( D, p6, c );
                    VecDot( p5c, c, &t6 );
                    //      + <\psi^(6) | D | \psi^(5)>
                    MatMult( D, p5, c );
                    VecDot( p6c, c, &t7 );
                    //      + <\psi^(7) | D | \psi^(4)>
                    MatMult( D, p4, c );
                    VecDot( p7c, c, &t8 );
                    //      + <\psi^(8) | D | \psi^(3)>
                    MatMult( D, p3, c );
                    VecDot( p8c, c, &t9 );
                    //      + <\psi^(9) | D | \psi^(2)>
                    MatMult( D, p2, c );
                    VecDot( p9c, c, &t10 );
                    //      + <\psi^(10) | D | \psi^(1)>
                    MatMult( D, p1, c );
                    VecDot( p10c, c, &t11 );
                    //      + <\psi^(11) | D | \psi^(0)>
                    MatMult( D, p0, c );
                    VecDot( p11c, c, &t12 );
                    if ( params.rank() == 0 )
                        std::cout << "terms: " << t1 << ", " << t2 << ", " << t3
                                  << ", " << t4 << ", " << t5 << ", " << t6
                                  << ", " << t7 << ", " << t8 << ", " << t9
                                  << ", " << t10 << ", " << t11 << ", " << t12
                                  << std::endl;
                    result += ( t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9 +
                                t10 + t11 + t12 );
                    VecDestroy( &p11 );
                    VecDestroy( &p10 );
                    VecDestroy( &p9 );
                    VecDestroy( &p8 );
                    VecDestroy( &p7 );
                    VecDestroy( &p6 );
                    VecDestroy( &p5 );
                    VecDestroy( &p4 );
                    VecDestroy( &p3 );
                    VecDestroy( &p2 );
                    VecDestroy( &p1 );
                    VecDestroy( &p0 );
                    VecDestroy( &p11c );
                    VecDestroy( &p10c );
                    VecDestroy( &p9c );
                    VecDestroy( &p8c );
                    VecDestroy( &p7c );
                    VecDestroy( &p6c );
                    VecDestroy( &p5c );
                    VecDestroy( &p4c );
                    VecDestroy( &p3c );
                    VecDestroy( &p2c );
                    VecDestroy( &p1c );
                    VecDestroy( &p0c );
                } while (
                    std::next_permutation( ( *i ).begin(), ( *i ).end() ) );

                if ( params.rank() == 0 )
                    std::cout << "final: " << result << std::endl;
                chi11_data[i - nparams.chi11s().begin()].push_back(
                    result * static_cast<double>( multiplicity ) /
                    static_cast<double>( math::factorial( 11 ) ) );
            }

            VecDestroy( &c );
        }
        H0 = params.read_energy_eigenvalues();
        if ( params.rank() == 0 ) std::cout << std::endl;
    }

    // export vectors.
    std::stringstream ss;
    for ( auto a = chi1_data.begin(); a != chi1_data.end(); ++a ) {
        ss << "chi1_" << nparams.chi1s()[a - chi1_data.begin()] << ".dat";
        common::export_vector_binary( ss.str(), *a );
        ss.str( "" );
    }
    ss.str( "" );

    for ( auto a = chi3_data.begin(); a != chi3_data.end(); ++a ) {
        ss << "chi3_";
        for ( auto b : nparams.chi3s()[a - chi3_data.begin()] ) ss << b << "_";
        ss << ".dat";
        common::export_vector_binary( ss.str(), *a );
        ss.str( "" );
    }

    ss.str( "" );
    for ( auto a = chi5_data.begin(); a != chi5_data.end(); ++a ) {
        ss << "chi5_";
        for ( auto b : nparams.chi5s()[a - chi5_data.begin()] ) ss << b << "_";
        ss << ".dat";
        common::export_vector_binary( ss.str(), *a );
        ss.str( "" );
    }
    ss.str( "" );
    for ( auto a = chi7_data.begin(); a != chi7_data.end(); ++a ) {
        ss << "chi7_";
        for ( auto b : nparams.chi7s()[a - chi7_data.begin()] ) ss << b << "_";
        ss << ".dat";
        common::export_vector_binary( ss.str(), *a );
        ss.str( "" );
    }
    ss.str( "" );
    for ( auto a = chi9_data.begin(); a != chi9_data.end(); ++a ) {
        ss << "chi9_";
        for ( auto b : nparams.chi9s()[a - chi9_data.begin()] ) ss << b << "_";
        ss << ".dat";
        common::export_vector_binary( ss.str(), *a );
        ss.str( "" );
    }
    ss.str( "" );
    for ( auto a = chi11_data.begin(); a != chi11_data.end(); ++a ) {
        ss << "chi11_";
        for ( auto b : nparams.chi11s()[a - chi11_data.begin()] )
            ss << b << "_";
        ss << ".dat";
        common::export_vector_binary( ss.str(), *a );
        ss.str( "" );
    }
    ss.str( "" );

    if ( params.rank() == 0 )
        std::cout << std::endl << "done with general code." << std::endl;


    PetscFinalize();
}

// Perturbative calculations:  \psi^(n) / E(w) etc.
// Each version takes "n" frequencies in a vector of real values

// template< typename Iterator >
Vec psi( int order, std::vector<double>::const_iterator frequencies_begin,
         std::vector<double>::const_iterator frequencies_end, PetscScalar wg,
         Vec& H0, Mat& D, Vec& psi0, Vec& mask,
         std::vector<BasisID>& /*prototype*/ )
{
    // we must have enough frequencies to do the calculation:
    assert( frequencies_end - frequencies_begin >= order );

    // the out vector (the one we return)
    Vec out;

    // a temporary vector
    Vec tmp;
    Vec orthog_projector;

    MPI_Comm comm;
    PetscObjectGetComm( (PetscObject)psi0, &comm );
    int rank;
    MPI_Comm_rank( comm, &rank );
    // make sure that out and tmp have the correct memory layout.
    VecDuplicate( psi0, &out );
    VecDuplicate( H0, &tmp );
    VecDuplicate( mask, orthog_projector );
    VecSet( orthog_projector, 1 );
    VecAXPY( orthog_projector, -1, mask );

    // out starts as psi0
    VecCopy( psi0, out );
    // VecCopy(H0, tmp);
    // VecSet(out, 1.);

    // do this "order" times
    for ( int i = 1; i <= order; ++i ) {
        // D \psi0 = tmp
        MatMult( D, out, tmp );
        // \psi0 = tmp
        VecCopy( tmp, out );
        // tmp = H0
        VecCopy( H0, tmp );
        // tmp = tmp - wg
        VecShift( tmp, -wg );

        // subtract the laser frequencies from H0 - wg
        for ( auto a = frequencies_begin + order - i;
              a < frequencies_begin + order; ++a ) {
            // tmp = tmp - \sum_a \omega_a
            VecShift( tmp, -( *a ) );
        }
        // tmp = 1/tmp
        VecReciprocal( tmp );

        // (1/tmp) * D | psi >
        VecPointwiseMult( out, tmp, out );
        VecPointwiseMult( out, orthog_projector, out );
    }

    // destroy the temporary
    VecDestroy( &tmp );

    return out;
}

Vec psi_conjugate( int                                 order,
                   std::vector<double>::const_iterator frequencies_begin,
                   std::vector<double>::const_iterator frequencies_end,
                   PetscScalar wg, Vec& H0, Mat& D, Vec& psi0, Vec& mask,
                   std::vector<BasisID>& /*prototype*/ )
{
    assert( frequencies_end - frequencies_begin >= order );
    Vec      out;
    Vec      tmp;
    Vec      orthog_projector;
    MPI_Comm comm;
    PetscObjectGetComm( (PetscObject)psi0, &comm );
    int rank;
    MPI_Comm_rank( comm, &rank );
    VecDuplicate( psi0, &out );
    VecDuplicate( H0, &tmp );
    VecCopy( psi0, out );

    VecDuplicate( mask, orthog_projector );
    VecSet( orthog_projector, 1 );
    VecAXPY( orthog_projector, -1, mask );

    for ( int i = 1; i <= order; ++i ) {
        MatMult( D, out, tmp );
        VecCopy( tmp, out );
        VecCopy( H0, tmp );
        VecShift( tmp, -wg );
        for ( auto a = frequencies_begin + i - 1; a >= frequencies_begin;
              --a ) {
            VecShift( tmp, ( *a ) );
        }
        VecReciprocal( tmp );
        VecPointwiseMult( out, tmp, out );

        VecPointwiseMult( out, orthog_projector, out );
    }

    VecDestroy( &tmp );

    return out;
}
