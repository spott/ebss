// mine:
#include <common/common.hpp>
#include <common/math.hpp>
#include <common/parameters/MomentumParameters.hpp>
#include <common/types.hpp>

// stl::
#include <vector>

// petsc:
#include <petsc.h>

// gsl:
#include <gsl/gsl_sf_legendre.h>


int main( int argc, const char** argv )
{
    int    ac = argc;
    char** av = new char*[argc];
    for ( size_t i = 0; i < argc; i++ ) {
        av[i] = new char[strlen( argv[i] ) + 1];
        std::copy( argv[i], argv[i] + strlen( argv[i] ) + 1, av[i] );
    }
    PetscInitialize( &ac, &av, PETSC_NULL, PETSC_NULL );

    MomentumParameters<double> kparams( argc, argv, MPI_COMM_WORLD );

    std::cout << kparams.print();

    if ( kparams.size() > 1 )
        std::cerr << "This program is incredibly inefficient / wrong if run on "
                     "more than one core"
                  << std::endl;

    Mat m = common::petsc_binary_read<Mat>( kparams.matrix_filename(),
                                            kparams.comm() );
    auto pro = kparams.kprototype();

    std::cout << "momentum prototype size: " << pro.size() << std::endl;
    int rows, cols;
    MatGetSize( m, &rows, &cols );
    std::cout << "matrix rows: " << rows << ", cols: " << cols << std::endl;
    auto sph_harmonics = math::gsl_spherical_harmonics(
        kparams.lmax() + 1, int( math::PI / kparams.dtheta() ) + 1 );

    Vec wf, k;
    for ( const auto& wf_fname : kparams.wf_filenames() ) {
        wf = common::petsc_binary_read<Vec>( wf_fname, kparams.comm() );
        VecGetSize( wf, &rows );
        std::cout << "wf " << wf_fname << " size: " << rows << std::endl;
        MatGetVecs( m, PETSC_NULL, &k );
        MatMult( m, wf, k );

        Matrix<std::complex<double>> spectrum(
            int( kparams.kmax() / kparams.dk() ),
            int( math::PI / kparams.dtheta() ) + 1 );
        std::cout << "spectrum: " << int( kparams.kmax() / kparams.dk() ) << ","
                  << int( math::PI / kparams.dtheta() ) + 1 << std::endl;

        // P (k, theta) = 2/pi | 1/k \sum_l i^l e^{i \sigma_l} Y_l (\theta)
        // \sum_n c_nl * <u_nl | u_kl> |^2
        // dP / dk_\rho dk_z = 2 pi k_rho P(k, \theta)

        // iterate over "l" for a specific k:
        PetscScalar* kval;
        VecGetArray( k, &kval );
        for ( size_t i = 0; i < pro.size(); ++i ) {
            for ( int t = 0;
                  t <= int( std::round( math::PI / kparams.dtheta() ) + .1 );
                  ++t ) {
                if ( pro[i].k > kparams.kmax() || pro[i].l > kparams.lmax() )
                    continue;
                spectrum( int( std::round( pro[i].k / kparams.dk() + .1 ) ) - 1,
                          t ) +=
                    std::sqrt( 2. / math::PI ) * 1. / pro[i].k *
                    std::pow( std::complex<double>( 0, 1 ), pro[i].l ) *
                    std::exp( std::complex<double>(
                        0, std::arg( math::Gamma_Lanczos( std::complex<double>(
                               pro[i].l + 1, -1. / pro[i].k ) ) ) ) ) *
                    sph_harmonics( t, pro[i].l ) * kval[i];
                // std::cout << int(pro[i].k/ kparams.dk()) << "," <<
                // kparams.dtheta() * t << ": "<<
                // spectrum(int(pro[i].k/kparams.dk()), t)  <<std::endl;
            }
        }

        Matrix<double> ddspectrum( spectrum.dim_x(), spectrum.dim_y() );
        for ( size_t x = 0; x < spectrum.dim_x(); ++x ) {
            for ( size_t y = 0; y < spectrum.dim_y(); ++y ) {
                ddspectrum( x, y ) =
                    std::pow( std::abs( spectrum( x, y ) ), 2 ) * 2 * math::PI *
                    ( x * kparams.dk() ) * std::sin( y * kparams.dtheta() );
            }
        }

        // write out!
        std::string f = wf_fname.substr( 0, wf_fname.size() - 4 );
        f.append( kparams.output_filename_postfix() );
        std::vector<int> prefix( {spectrum.dim_x(), spectrum.dim_y()} );
        // prefix[0] = * (reinterpret_cast<char*>(new int[2]{spectrum.dim_x(),
        // spectrum.dim_y()} ));
        // prefix.data() +  = reinterpret_cast<char>( spectrum.dim_x() );
        common::export_vector_binary( f, ddspectrum.data(), prefix );
    }
    VecDestroy( &wf );
    VecDestroy( &k );
    MatDestroy( &m );
    PetscFinalize();
}
