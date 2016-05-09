#pragma once

#include <array>
#include <csignal>
#include <cassert>
#include <slepceps.h>

typedef struct {
    Mat* D;
    Vec* H;
    HamiltonianParameters<PetscReal>* hparams;
    PulsetrainParameters* laser;
    AbsorberParameters* absorber;
    DipoleParameters* dipole;
} context;

namespace cranknicholson
{

volatile sig_atomic_t sig = 0;

// PetscErrorCode signal_handler( int signal, void*) {
//     std::cerr << "Caught signal: " << signal << "... " << std::endl;
//     sig = signal; }
void signal_handler( int signal) {
    std::cerr << "Caught signal: " << signal << "... " << std::endl;
    sig = signal; }

PetscErrorCode solve( Vec* wf, context* cntx, Mat* A)
{
    // PetscPushSignalHandler(&signal_handler, NULL);
    bool restore = true;
    PetscInt zero = 0;
    try {
        std::cout << "looking for recovery file" << std::endl;
        zero = common::import_vector_binary<int>("./failed").back();
        std::cout << "found recovery file, got: " << zero << " from it" << std::endl;
    } catch (...) {
        std::cout << "didn't find a recovery file.  proceeding as if from scratch" << std::endl;
        restore = false;
        zero = 0;
    }
    PetscViewer view;
    KSP ksp;
    PC pc;
    EPS eps;
    EPSCreate( cntx->hparams->comm(), &eps );
    KSPCreate( cntx->hparams->comm(), &ksp );
    KSPSetType( ksp, KSPGMRES );
    KSPGetPC( ksp, &pc );
    PCSetType( pc, PCJACOBI );
    KSPSetTolerances( ksp, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT );
    KSPSetFromOptions( ksp );

    std::vector<PetscReal> efvec;
    std::vector<PetscReal> time;

    PetscScalar cn_factor = std::complex<double>( 0, -0.5 * cntx->laser->dt() );
    PetscReal t = 0;
    PetscScalar ef = cntx->laser->efield( t );
    PetscReal maxtime = cntx->laser->max_time();
    PetscInt step = 0;

    if (restore)
        {
            if ( cntx->hparams->rank() == 0 ) {
                std::cout << "getting time and ef " << std::endl;
                time = common::import_vector_binary<PetscReal>( "time.dat");
                efvec = common::import_vector_binary<PetscReal>(cntx->laser->laser_filename());
                t = time.back();
                ef = efvec.back();
                step = time.size();
            }
            else
                {
                    std::cout << "getting time and ef " << std::endl;
                    auto temp_time = common::import_vector_binary<PetscReal>( "time.dat");
                    auto temp_efvec = common::import_vector_binary<PetscReal>(cntx->laser->laser_filename());
                    t = temp_time.back();
                    ef = temp_efvec.back();
                    step = temp_time.size();
                }

            std::cout << "getting wf " << std::endl;
            std::string filename = std::string("wf_" + std::to_string(zero) + ".dat");
            VecDestroy(wf);
            *wf = common::petsc_binary_read<Vec>(filename, cntx->hparams->comm());

        }

    Vec tmp;
    Vec prob;
    PetscReal norm;
    // PetscInt        p = 0;
    // PetscReal       period = maxtime / (8 * cntx->laser->cycles()) ;
    // PetscReal       p0 = period;

    if ( cntx->hparams->rank() == 0 )
        std::cout << "maxtime: " << maxtime
                  << " steps: " << maxtime / cntx->laser->dt()
                  << " pulse length: " << cntx->laser->pulse_length()
                  << std::endl;

    VecDuplicate( *wf, &tmp );
    VecAssemblyBegin( tmp );
    VecAssemblyEnd( tmp );

    Vec abs;
    VecDuplicate( *wf, &abs );
    cntx->absorber->absorb( &abs, cntx->hparams );
    // if (cntx->absorber->type() == "cx_rot" || cntx->absorber->type() ==
    // "cx_scale")
    // VecPointwiseMult(*(cntx->H), *(cntx->H), abs);

    std::string file_name = std::string( "./absorber.dat" );
    PetscViewerASCIIOpen( cntx->hparams->comm(), file_name.c_str(), &view );
    PetscViewerSetFormat( view, PETSC_VIEWER_ASCII_SYMMODU );
    VecView( abs, view );

    VecDuplicate( *wf, &prob );
    VecAssemblyBegin( prob );
    VecAssemblyEnd( prob );

    std::vector<std::ofstream*> dipole;
    std::ofstream population( "instant_pop.csv", restore ? std::ios::ate : std::ios::out );
    // dipol( 3 * cntx->dipole->decompositions().size() + 1 );
    if ( cntx->hparams->rank() == 0 ) {
        assert( cntx->dipole->dipole_filename().size() > 0 );
        for ( auto& a : cntx->dipole->dipole_filename() ) {
            try {
                if (restore)
                    {
                        dipole.emplace_back(
                            new std::ofstream( a, std::ios::binary | std::ios::ate ) );
                        dipole.back()->seekp( sizeof(PetscScalar) * step );
                    }
                else
                    dipole.emplace_back(
                        new std::ofstream( a, std::ios::binary | std::ios::out ) );
            } catch ( ... ) {
                std::cerr << "couldn't open the dipole file, oops... it "
                             "won't be "
                             "written to disk" << std::endl;
            }
        }
        // population = std::move(std::ofstream( "instant_pop.csv",
        // std::ios::ate ));
    }
    // std::vector< std::vector<PetscReal> > dipole( 3 *
    // cntx->dipole->decompositions().size() + 1 );

    PetscReal norm_lost = 0;

    std::signal( SIGUSR1, signal_handler );

    Vec eigenvector;
    MatGetVecs( *( cntx->D ), PETSC_NULL, &eigenvector );

    PetscScalar eff = ef;
    // bool at_zero = false;
    double dt_ = cntx->laser->dt();
    while ( t < maxtime ) {

        auto ef = eff;
        // ef-forward:
        t += cntx->laser->dt();
        eff = cntx->laser->efield( t );
        //		if (step > 2 && math::signum(eff.real()) !=
        //math::signum(ef.real()))
        //		{
        //			auto ddt = dt_ / 2.;
        //			while( std::abs(eff ) > 1e-15 and (math::signum(eff.real()) ==
        //math::signum(ef.real())) )
        //			{
        //				if (math::signum(eff.real()) != math::signum(ef.real()))
        //				{
        //					dt_ -= ddt;
        //				}
        //				else if (math::signum(eff.real()) ==
        //math::signum(ef.real()))
        //				{
        //					dt_ += ddt;
        //				}
        //				eff = cntx->laser->efield(t + dt_);
        //				ddt /= 2.;
        //			}
        //			t = t + dt_;
        //			cn_factor = std::complex<double>(0, -.5 * dt_);
        //			dt_ = cntx->laser->dt();
        //		}
        //		else
        //		{
        //			t += dt_;
        //			cn_factor = std::complex<double>(0, -.5 * dt_);
        //		}
        MatCopy( *( cntx->D ), *A, SAME_NONZERO_PATTERN ); // A = D
        // This has different 't's on both sides:
        MatScale( *A, ( ef ) ); // A = ef(t) * D
        MatDiagonalSet( *A, *( cntx->H ),
                        INSERT_VALUES ); // A = ef(t) * D + H_0

        MatScale( *A, cn_factor ); // A = - i * .5 * dt [ ef(t) * D + H_0 ]
        MatShift( *A, std::complex<double>( 1, 0 ) ); // A = - i * .5 * dt [
                                                      // ef(t) * D + H_0 ] + 1
        MatMult( *A, *wf, tmp );                      // A u_n = tmp
        MatCopy( *( cntx->D ), *A, SAME_NONZERO_PATTERN ); // A = D
        // This has different 't's on both sides:
        MatScale( *A, ( eff ) ); // A = eff(t) * D
        MatDiagonalSet( *A, *( cntx->H ),
                        INSERT_VALUES ); // A = eff(t) * D + H_0
        MatScale( *A, -cn_factor );      // A = i * .5 dt[ef(t+dt) * D + H_0 ]
        MatShift( *A, std::complex<double>( 1, 0 ) ); // A = i * .5 dt [ef(t+dt)
                                                      // * D + H_0 ] + 1


        KSPSetOperators( ksp, *A, *A ); //, SAME_NONZERO_PATTERN); // Solve[ A x
                                        // = tmp ]
                                        // for x
        KSPSolve( ksp, tmp, *wf );

        if ( cntx->absorber->type() == "cosine" ) {

            PetscReal n = 0;
            PetscReal n2 = 0;
            VecNorm( *wf, NORM_2, &n );
            VecPointwiseMult( *wf, *wf, abs );
            VecNorm( *wf, NORM_2, &n2 );
            norm_lost += n2 - n;
            // if (cntx->hparams->rank() == 0) std::cout << "lost " << n2-n
            // << "
            // norm, total = " << norm_lost << std::endl;
        }

        // look at the next point
        // if we would be at a zero, search for the ACTUAL zero:

        // check if we are "in" a pulse now:
        if ( !cntx->laser->in_pulse( t ) ) {
            // first write out the zero:
            std::cout << output::red;
            if ( cntx->hparams->rank() == 0 )
                std::cout << "time: " << t << " step: " << step
                          << " efield: " << ef << " norm-1: " << norm - 1
                          << " * " << zero << std::endl;
            std::ostringstream wf_name;
            wf_name << "./wf_" << zero << ".dat";
            PetscViewerBinaryOpen( cntx->hparams->comm(), wf_name.str().c_str(),
                                   FILE_MODE_WRITE, &view );
            zero++;
            // analytically propagate for said time:
            std::cout << "analytically propagating the pulse: ";
            PetscReal dt = cntx->laser->next_pulse_start( t ) -
                           t; // the difference between then and now.
            std::cout << dt << std::endl;
            math::FieldFreePropagate( ( cntx->H ), wf,
                                      dt ); // propagate forward in time
            t = cntx->laser->next_pulse_start( t );
            ef = cntx->laser->efield( t );
            // write out the next zero:
            if ( cntx->hparams->rank() == 0 )
                std::cout << "time: " << t << " step: " << step
                          << " efield: " << eff.real()
                          << " norm-1: " << norm - 1 << " *" << zero
                          << std::endl;
            wf_name.str( "" );
            wf_name << "./wf_" << zero << ".dat";
            PetscViewerBinaryOpen( cntx->hparams->comm(), wf_name.str().c_str(),
                                   FILE_MODE_WRITE, &view );
            // VecView( *wf, view );
            zero++;
            std::cout << output::reset;
            continue;
        }
        // else do nothing...

        step++;
        // at the zero of the field: write out the vector:
        if ( math::signum( ef.real() ) != math::signum( eff.real() ) ) {
            if ( cntx->hparams->rank() == 0 )
                {
                std::cout << "time: " << t << " step: " << step
                          << " efield: " << eff.real()
                          << " norm-1: " << norm - 1 << " *" << std::endl;
                if (cntx->hparams->rank() == 0)
                    {
                        auto tvec = std::vector<int>();
                        tvec.push_back(zero);
                        common::export_vector_binary("failed",tvec);
                    }
                if ( cntx->hparams->rank() == 0 ) {
                    try {
                        common::export_vector_binary( "time.dat", time );
                    } catch ( ... ) {
                        std::cerr << "couldn't open the time file, oops... it won't be "
                            "written to disk" << std::endl;
                    }
                    try {
                        common::export_vector_binary( cntx->laser->laser_filename(),
                                                      efvec );
                    } catch ( ... ) {
                        std::cerr << "couldn't open the efvec file, oops... it won't be "
                            "written to disk" << std::endl;
                    }
                    for (auto& a: dipole)
                        a->flush();
                }
                }
            std::ostringstream wf_name;
            wf_name << "./wf_" << zero << ".dat";
            PetscViewerBinaryOpen( cntx->hparams->comm(), wf_name.str().c_str(),
                                   FILE_MODE_WRITE, &view );
            VecView( *wf, view );
            zero++;
        }

        if ( true ) //(!(step%10))
        {
            if ( cntx->hparams->rank() == 0 ) {
                time.push_back( t );
                efvec.push_back( eff.real() );
            }
            cntx->dipole->find_dipole_moment_decompositions(
                *( cntx->D ), *wf, dipole, cntx->hparams->prototype(),
                ( std::pow( cntx->laser->envelope( t, 0 ), 2 ) *
                  cntx->laser->intensity() /
                  ( 4. * cntx->laser->frequency() * cntx->laser->frequency() ) )
                    .real() );
        }
        if ( !( step % 100 ) || step < 100 ) {
            VecCopy( *wf, prob );
            VecAbs( prob );
            VecNorm( prob, NORM_2, &norm );
            if ( cntx->hparams->rank() == 0 )
                std::cout << "time: " << t << " step: " << step
                          << " efield: " << eff.real()
                          << " norm-1: " << norm - 1
                          << " norm lost to absorbers: " << norm_lost
                          << std::endl;
        } else if ( cntx->hparams->rank() == 0 )
            std::cout << "time: " << t << " step: " << step
                      << " efield: " << eff.real() << std::endl;


        if ( sig > 0 )
            {
                //make recovery vector, write out "zero number"
                std::cerr << "got interrupt signal " << sig << " breaking..." << std::endl;
                if (cntx->hparams->rank() == 0)
                    {
                        auto tvec = std::vector<int>();
                        tvec.push_back(zero);
                        common::export_vector_binary("failed",tvec);
                    }
                break;
                //write out wf:
                //write out time vector:
                //write out efield vector:
                // close dipole files
            }
    }
    //MPI_Barrier(cntx->hparams->comm());
    file_name = std::string( "./wf_final.dat" );
    if (sig>0)
        file_name = std::string( "./wf_interupted.dat");
    PetscViewerBinaryOpen( cntx->hparams->comm(), file_name.c_str(),
                           FILE_MODE_WRITE, &view );
    VecView( *wf, view );

    // we don't want to crash the nodes from lack of memory, so we only do
    // this on the head node
    if ( cntx->hparams->rank() == 0 ) {
        try {
            common::export_vector_binary( "time.dat", time );
        } catch ( ... ) {
            std::cerr << "couldn't open the time file, oops... it won't be "
                         "written to disk" << std::endl;
        }
        try {
            common::export_vector_binary( cntx->laser->laser_filename(),
                                          efvec );
        } catch ( ... ) {
            std::cerr << "couldn't open the efvec file, oops... it won't be "
                         "written to disk" << std::endl;
        }

        // close the dipole files
        for ( auto& a : dipole ) {
            a->close();
        }
    }

    KSPDestroy( &ksp );
    std::cerr << "leaving cranknicholson" << std::endl;
    return 0;
}
}
