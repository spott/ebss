
#include <petsc.h>

// stl:
#include <complex>
#include <vector>
#include <cassert>
#include <utility>

// mine:
#include <common/parameters/HamiltonianParameters.hpp>
#include <common/parameters/NonlinearParameters.hpp>
#include <common/parameters/StateParameters.hpp>
#include <common/math.hpp>
#include <common/common.hpp>


Vec psi( int order,
         std::vector<double>::const_iterator frequencies_begin,
         std::vector<double>::const_iterator frequencies_end,
         PetscScalar wg,
         Vec& H0,
         Mat& D,
         Vec& psi0,
         Vec& mask,
         std::vector<BasisID>& prototype );
Vec psi_conjugate( int order,
                   std::vector<double>::const_iterator frequencies_begin,
                   std::vector<double>::const_iterator frequencies_end,
                   PetscScalar wg,
                   Vec& H0,
                   Mat& D,
                   Vec& psi0,
                   Vec& mask,
                   std::vector<BasisID>& prototype );

struct VectorElement
{
    const std::tuple<std::complex<double>, int>& el;
    const std::vector<BasisID>& proto;

    VectorElement( const std::tuple<std::complex<double>, int>& element,
                   const std::vector<BasisID>& prototype )
        : el( element ), proto( prototype ) {};
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
    int ac = argc;
    char** av = new char* [argc];
    for ( size_t i = 0; i < argc; i++ ) {
        av[i] = new char[strlen( argv[i] ) + 1];
        std::copy( argv[i], argv[i] + strlen( argv[i] ) + 1, av[i] );
    }
    PetscInitialize( &ac, &av, PETSC_NULL, PETSC_NULL );

    // the parameters (where to find the hamiltonian)
    NonlinearParameters nparams( argc, argv, MPI_COMM_WORLD );
    HamiltonianParameters<double> params( argc, argv, MPI_COMM_WORLD );
    StateParameters sparams( argc, argv, MPI_COMM_WORLD );

    auto prototype = params.prototype();

    // the parameters for this:
    // std::cout << nparams.print();

    std::cout << std::scientific;

    // read in the matrices
    Mat D = params.read_dipole_matrix();
    Vec H0 = params.read_energy_eigenvalues();
    std::cout << prototype.size() << std::endl;
    int size;
    VecGetSize( H0, &size );
    std::cout << size << std::endl;

    assert( size == prototype.size() );
    std::function<bool(int, int)> dipole_selection_rules = [prototype](
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
    Vec mask = common::map_function( H0, []( PetscScalar a, PetscInt in ) {
        return ( std::abs( a.real() + .5 ) <= 1e-4 ) ? 0 : 1;
    } );
    // Vec mask = common::map_function(H0, [](PetscScalar in) { return 1; });

    auto imgs = nparams.imgs();
    // VecShift(H0, std::complex<double>(0,-.00001));
    Vec psi0;
    VecDuplicate( H0, &psi0 );
    VecSetValue( psi0, 0, 1., INSERT_VALUES );
    VecAssemblyBegin( psi0 );
    VecAssemblyEnd( psi0 );

    // find the resonances:
    PetscScalar wg;
    VecDot( H0, psi0, &wg );
    Vec resonances;
    VecDuplicate( H0, &resonances );
    VecSet( resonances, wg );
    VecAXPY( resonances, -1., H0 ); // single resonances

    PetscViewer a;
    PetscViewerASCIIOpen( params.comm(), "resonances", &a );
    PetscViewerSetFormat( a, PETSC_VIEWER_ASCII_INDEX );

    // get the frequencies list from nonlinear params:
    auto freqs = nparams.freqs();
    std::sort(freqs.begin(), freqs.end());

    try
    {
        common::export_vector_binary( std::string( "./frequencies.dat" ),
                                      freqs );
    }
    catch ( ... )
    {
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

        for ( int f = 0; f < freqs.size(); ++f ) {
            if ( params.rank() == 0 )
                std::cout << f << "(" << freqs[f] << ")" << std::endl;
            PetscScalar t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
            Vec c;
            VecDuplicate( H0, &c );
            // we need to loop through all permutations of all chi types:

            // first we need to generate the correct vectors:
            // Chi1
            for ( auto i = nparams.chi1s().begin(); i != nparams.chi1s().end();
                  ++i ) {
                if ( params.rank() == 0 )
                    std::cout << "========================================="
                              << std::endl << *i << std::endl;
                std::vector<double> freq{( *i ) * freqs[f]};
                Vec p0 = psi( 0,
                              freq.cbegin(),
                              freq.cend(),
                              wg,
                              H0,
                              D,
                              psi0,
                              mask,
                              prototype );
                Vec p1 = psi( 1,
                              freq.cbegin(),
                              freq.cend(),
                              wg,
                              H0,
                              D,
                              psi0,
                              mask,
                              prototype );
                Vec p1c = psi_conjugate( 1,
                                         freq.cbegin(),
                                         freq.cend(),
                                         wg,
                                         H0,
                                         D,
                                         psi0,
                                         mask,
                                         prototype );
                Vec p0c = psi_conjugate( 0,
                                         freq.cbegin(),
                                         freq.cend(),
                                         wg,
                                         H0,
                                         D,
                                         psi0,
                                         mask,
                                         prototype );

                // chi1 = <\psi^(0) | D | \psi^(1)>
                MatMult( D, p1, c );
                VecDot( p0c, c, &t1 );
                // chi1 = <\psi^(1) | D | \psi^(0)>
                MatMult( D, p0, c );
                VecDot( p1c, c, &t2 );
                chi1_data[i - nparams.chi1s().begin()].push_back( ( t1 + t2 ) );
                if ( params.rank() == 0 )
                    std::cout << "terms: " << t1 << ", " << t2 << std::endl;
                if ( params.rank() == 0 )
                    std::cout << "final: " << ( t1 + t2 ) << std::endl;
                VecDestroy( &p1 );
                VecDestroy( &p0 );
                VecDestroy( &p1c );
                VecDestroy( &p0c );
            }

            // Chi3
            for ( auto i = nparams.chi3s().begin(); i != nparams.chi3s().end();
                  ++i ) {
                if ( params.rank() == 0 )
                    std::cout << "========================================="
                              << std::endl << *i << std::endl;
                std::sort( ( *i ).begin(), ( *i ).end() );
                size_t multiplicity = 1;
                std::array<int, 3> ts{{0, 0, 0}};
                for ( auto m : ( *i ) ) {
                    if ( m == -1 ) ts[0]++;
                    if ( m == 1 ) ts[2]++;
                    if ( m == 0 ) ts[1]++;
                }
                for ( auto m : ts ) multiplicity *= math::factorial( m );

                PetscScalar result;
                do {

                    std::vector<double> freq{( *i )[0] * freqs[f],
                                             ( *i )[1] * freqs[f],
                                             ( *i )[2] * freqs[f]};
                    Vec p3 = psi( 3,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p2 = psi( 2,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p0 = psi( 0,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p1 = psi( 1,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p3c = psi_conjugate( 3,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p2c = psi_conjugate( 2,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p1c = psi_conjugate( 1,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p0c = psi_conjugate( 0,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );

                    // chi3 = <\psi^(0) | D | \psi^(3)>
                    MatMult( D, p3, c );
                    VecDot( p0c, c, &t1 );
                    //      + <\psi^(1) | D | \psi^(2)>
                    MatMult( D, p2, c );
                    VecDot( p1c, c, &t2 );
                    //      + <\psi^(2) | D | \psi^(1)>
                    MatMult( D, p1, c );
                    VecDot( p2c, c, &t3 );
                    ////      + <\psi^(3) | D | \psi^(0)>
                    MatMult( D, p0, c );
                    VecDot( p3c, c, &t4 );
                    result += ( t1 + t2 + t3 + t4 );
                    if ( params.rank() == 0 )
                        std::cout << "terms: " << t1 << ", " << t2 << ", " << t3
                                  << ", " << t4 << std::endl;
                    VecDestroy( &p3 );
                    VecDestroy( &p2 );
                    VecDestroy( &p1 );
                    VecDestroy( &p0 );
                    VecDestroy( &p3c );
                    VecDestroy( &p2c );
                    VecDestroy( &p1c );
                    VecDestroy( &p0c );
                } while (
                    std::next_permutation( ( *i ).begin(), ( *i ).end() ) );
                if ( params.rank() == 0 )
                    std::cout << "final: " << result << std::endl;
                chi3_data[i - nparams.chi3s().begin()]
                    .push_back( result * static_cast<double>( multiplicity ) /
                                static_cast<double>( math::factorial( 3 ) ) );
            }
            // chi5
            for ( auto i = nparams.chi5s().begin(); i != nparams.chi5s().end();
                  ++i ) {
                if ( params.rank() == 0 )
                    std::cout << "========================================="
                              << std::endl << *i << std::endl;
                std::sort( ( *i ).begin(), ( *i ).end() );
                size_t multiplicity = 1;
                std::array<int, 3> ts{ {0, 0, 0} };
                for ( auto m : ( *i ) ) {
                    if ( m == -1 ) ts[0]++;
                    if ( m == 0 ) ts[1]++;
                    if ( m == 1 ) ts[2]++;
                }
                for ( auto m : ts ) multiplicity *= math::factorial( m );
                PetscScalar result;
                do {
                    std::vector<double> freq{
                        ( *i )[0] * freqs[f], ( *i )[1] * freqs[f],
                        ( *i )[2] * freqs[f], ( *i )[3] * freqs[f],
                        ( *i )[4] * freqs[f]};
                    // if (params.rank() == 0) std::cout << "p5: " << std::endl;
                    Vec p5 = psi( 5,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    // if (params.rank() == 0) std::cout << "p4: " << std::endl;
                    Vec p4 = psi( 4,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    // if (params.rank() == 0) std::cout << "p3: " << std::endl;
                    Vec p3 = psi( 3,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    // if (params.rank() == 0) std::cout << "p2: " << std::endl;
                    Vec p2 = psi( 2,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    // if (params.rank() == 0) std::cout << "p1: " << std::endl;
                    Vec p1 = psi( 1,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    // if (params.rank() == 0) std::cout << "p0: " << std::endl;
                    Vec p0 = psi( 0,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    // if (params.rank() == 0) std::cout << "p5c: " <<
                    // std::endl;
                    Vec p5c = psi_conjugate( 5,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    // if (params.rank() == 0) std::cout << "p4c: " <<
                    // std::endl;
                    Vec p4c = psi_conjugate( 4,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    // if (params.rank() == 0) std::cout << "p3c: " <<
                    // std::endl;
                    Vec p3c = psi_conjugate( 3,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    // if (params.rank() == 0) std::cout << "p2c: " <<
                    // std::endl;
                    Vec p2c = psi_conjugate( 2,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    // if (params.rank() == 0) std::cout << "p1c: " <<
                    // std::endl;
                    Vec p1c = psi_conjugate( 1,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    // if (params.rank() == 0) std::cout << "p0c: " <<
                    // std::endl;
                    Vec p0c = psi_conjugate( 0,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );

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
                chi5_data[i - nparams.chi5s().begin()]
                    .push_back( result * static_cast<double>( multiplicity ) /
                                static_cast<double>( math::factorial( 5 ) ) );
            }
            // Chi7
            for ( auto i = nparams.chi7s().begin(); i != nparams.chi7s().end();
                  ++i ) {
                if ( params.rank() == 0 )
                    std::cout << "========================================="
                              << std::endl << *i << std::endl;
                std::sort( ( *i ).begin(), ( *i ).end() );
                size_t multiplicity = 1;
                std::array<int, 3> ts{ {0, 0, 0} };
                for ( auto m : ( *i ) ) {
                    if ( m == -1 ) ts[0]++;
                    if ( m == 0 ) ts[1]++;
                    if ( m == 1 ) ts[2]++;
                }
                for ( auto m : ts ) multiplicity *= math::factorial( m );
                PetscScalar result;
                do {
                    std::vector<double> freq{
                        ( *i )[0] * freqs[f], ( *i )[1] * freqs[f],
                        ( *i )[2] * freqs[f], ( *i )[3] * freqs[f],
                        ( *i )[4] * freqs[f], ( *i )[5] * freqs[f],
                        ( *i )[6] * freqs[f]};
                    Vec p7 = psi( 7,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p6 = psi( 6,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p5 = psi( 5,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p4 = psi( 4,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p3 = psi( 3,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p2 = psi( 2,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p1 = psi( 1,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p0 = psi( 0,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p7c = psi_conjugate( 7,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p6c = psi_conjugate( 6,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p5c = psi_conjugate( 5,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p4c = psi_conjugate( 4,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p3c = psi_conjugate( 3,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p2c = psi_conjugate( 2,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p1c = psi_conjugate( 1,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p0c = psi_conjugate( 0,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );

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
                chi7_data[i - nparams.chi7s().begin()]
                    .push_back( result * static_cast<double>( multiplicity ) /
                                static_cast<double>( math::factorial( 7 ) ) );
            }

            // Chi9
            for ( auto i = nparams.chi9s().begin(); i != nparams.chi9s().end();
                  ++i ) {
                if ( params.rank() == 0 )
                    std::cout << "========================================="
                              << std::endl << *i << std::endl;
                std::sort( ( *i ).begin(), ( *i ).end() );
                size_t multiplicity = 1;
                std::array<int, 3> ts{ {0, 0, 0} };
                for ( auto m : ( *i ) ) {
                    if ( m == -1 ) ts[0]++;
                    if ( m == 0 ) ts[1]++;
                    if ( m == 1 ) ts[2]++;
                }
                for ( auto m : ts ) multiplicity *= math::factorial( m );
                PetscScalar result;
                do {
                    std::vector<double> freq{
                        ( *i )[0] * freqs[f], ( *i )[1] * freqs[f],
                        ( *i )[2] * freqs[f], ( *i )[3] * freqs[f],
                        ( *i )[4] * freqs[f], ( *i )[5] * freqs[f],
                        ( *i )[6] * freqs[f], ( *i )[7] * freqs[f],
                        ( *i )[8] * freqs[f]};
                    Vec p9 = psi( 9,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p8 = psi( 8,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p7 = psi( 7,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p6 = psi( 6,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p5 = psi( 5,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p4 = psi( 4,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p3 = psi( 3,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p2 = psi( 2,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p1 = psi( 1,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p0 = psi( 0,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p9c = psi_conjugate( 9,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p8c = psi_conjugate( 8,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p7c = psi_conjugate( 7,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p6c = psi_conjugate( 6,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p5c = psi_conjugate( 5,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p4c = psi_conjugate( 4,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p3c = psi_conjugate( 3,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p2c = psi_conjugate( 2,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p1c = psi_conjugate( 1,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p0c = psi_conjugate( 0,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );

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
                chi9_data[i - nparams.chi9s().begin()]
                    .push_back( result * static_cast<double>( multiplicity ) /
                                static_cast<double>( math::factorial( 9 ) ) );
            }

            // Chi11
            for ( auto i = nparams.chi11s().begin();
                  i != nparams.chi11s().end();
                  ++i ) {
                if ( params.rank() == 0 )
                    std::cout << "========================================="
                              << std::endl << *i << std::endl;
                std::sort( ( *i ).begin(), ( *i ).end() );
                size_t multiplicity = 1;
                std::array<int, 3> ts{ {0, 0, 0} };
                for ( auto m : ( *i ) ) {
                    if ( m == -1 ) ts[0]++;
                    if ( m == 0 ) ts[1]++;
                    if ( m == 1 ) ts[2]++;
                }
                for ( auto m : ts ) multiplicity *= math::factorial( m );
                PetscScalar result;
                do {
                    std::vector<double> freq{
                        ( *i )[0] * freqs[f], ( *i )[1] * freqs[f],
                        ( *i )[2] * freqs[f], ( *i )[3] * freqs[f],
                        ( *i )[4] * freqs[f], ( *i )[5] * freqs[f],
                        ( *i )[6] * freqs[f], ( *i )[7] * freqs[f],
                        ( *i )[8] * freqs[f], ( *i )[9] * freqs[f],
                        ( *i )[10] * freqs[f]};

                    Vec p11 = psi( 11,
                                   freq.cbegin(),
                                   freq.cend(),
                                   wg,
                                   H0,
                                   D,
                                   psi0,
                                   mask,
                                   prototype );
                    Vec p10 = psi( 10,
                                   freq.cbegin(),
                                   freq.cend(),
                                   wg,
                                   H0,
                                   D,
                                   psi0,
                                   mask,
                                   prototype );
                    Vec p9 = psi( 9,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p8 = psi( 8,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p7 = psi( 7,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p6 = psi( 6,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p5 = psi( 5,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p4 = psi( 4,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p3 = psi( 3,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p2 = psi( 2,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p1 = psi( 1,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p0 = psi( 0,
                                  freq.cbegin(),
                                  freq.cend(),
                                  wg,
                                  H0,
                                  D,
                                  psi0,
                                  mask,
                                  prototype );
                    Vec p11c = psi_conjugate( 11,
                                              freq.cbegin(),
                                              freq.cend(),
                                              wg,
                                              H0,
                                              D,
                                              psi0,
                                              mask,
                                              prototype );
                    Vec p10c = psi_conjugate( 10,
                                              freq.cbegin(),
                                              freq.cend(),
                                              wg,
                                              H0,
                                              D,
                                              psi0,
                                              mask,
                                              prototype );
                    Vec p9c = psi_conjugate( 9,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p8c = psi_conjugate( 8,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p7c = psi_conjugate( 7,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p6c = psi_conjugate( 6,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p5c = psi_conjugate( 5,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p4c = psi_conjugate( 4,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p3c = psi_conjugate( 3,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p2c = psi_conjugate( 2,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p1c = psi_conjugate( 1,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );
                    Vec p0c = psi_conjugate( 0,
                                             freq.cbegin(),
                                             freq.cend(),
                                             wg,
                                             H0,
                                             D,
                                             psi0,
                                             mask,
                                             prototype );

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
                chi11_data[i - nparams.chi11s().begin()]
                    .push_back( result * static_cast<double>( multiplicity ) /
                                static_cast<double>( math::factorial( 11 ) ) );
            }

            VecDestroy( &c );
        }
        H0 = params.read_energy_eigenvalues();
        std::cout << std::endl;
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

    std::cout << std::endl << "done with general code." << std::endl;


    PetscFinalize();
}

// Perturbative calculations:  \psi^(n) / E(w) etc.
// Each version takes "n" frequencies in a vector of real values

// template< typename Iterator >
Vec psi( int order,
         std::vector<double>::const_iterator frequencies_begin,
         std::vector<double>::const_iterator frequencies_end,
         PetscScalar wg,
         Vec& H0,
         Mat& D,
         Vec& psi0,
         Vec& mask,
         std::vector<BasisID>& prototype )
{
    // we must have enough frequencies to do the calculation:
    assert( frequencies_end - frequencies_begin >= order );

    // the out vector (the one we return)
    Vec out;

    // a temporary vector
    Vec tmp;

    MPI_Comm comm;
    PetscObjectGetComm( (PetscObject)psi0, &comm );
    int rank;
    MPI_Comm_rank( comm, &rank );
    // make sure that out and tmp have the correct memory layout.
    VecDuplicate( psi0, &out );
    VecDuplicate( H0, &tmp );

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
              a < frequencies_begin + order;
              ++a ) {
            // tmp = tmp - \sum_a \omega_a
            VecShift( tmp, -( *a ) );
        }
        // tmp = 1/tmp
        VecReciprocal( tmp );
        VecPointwiseMult( out, tmp, out );
        VecPointwiseMult( out, mask, out );
    }

    // destroy the temporary
    VecDestroy( &tmp );

    return out;
}

Vec psi_conjugate( int order,
                   std::vector<double>::const_iterator frequencies_begin,
                   std::vector<double>::const_iterator frequencies_end,
                   PetscScalar wg,
                   Vec& H0,
                   Mat& D,
                   Vec& psi0,
                   Vec& mask,
                   std::vector<BasisID>& prototype )
{
    assert( frequencies_end - frequencies_begin >= order );
    Vec out;
    Vec tmp;
    MPI_Comm comm;
    PetscObjectGetComm( (PetscObject)psi0, &comm );
    int rank;
    MPI_Comm_rank( comm, &rank );
    VecDuplicate( psi0, &out );
    VecDuplicate( H0, &tmp );
    VecCopy( psi0, out );

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
        VecPointwiseMult( out, mask, out );

    }

    VecDestroy( &tmp );

    return out;
}
