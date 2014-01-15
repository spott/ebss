#pragma once

#include <array>

typedef struct
{
    Mat* D;
    Vec* H;
    HamiltonianParameters<PetscReal>* hparams;
    PulsetrainParameters* laser;
    AbsorberParameters* absorber;
    DipoleParameters* dipole;
} context;

namespace cranknicholson
{

PetscErrorCode solve( Vec* wf, context* cntx, Mat* A )
{
    PetscViewer view;
    KSP ksp;
    PC pc;
    KSPCreate( cntx->hparams->comm(), &ksp );
    KSPSetType( ksp, KSPGMRES );
    KSPGetPC( ksp, &pc );
    PCSetType( pc, PCJACOBI );
    KSPSetTolerances(
        ksp, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT );
    KSPSetFromOptions( ksp );

    PetscScalar cn_factor =
        std::complex<double>( 0, -0.5 * cntx->laser->dt() );
    PetscReal t = 0;
    PetscScalar ef = cntx->laser->efield( t );
    PetscReal maxtime = cntx->laser->max_time();
    PetscInt step = 0;
    Vec tmp;
    Vec prob;
    PetscReal norm;
    PetscInt zero = 0;
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
    PetscViewerASCIIOpen(
        cntx->hparams->comm(), file_name.c_str(), &view );
    PetscViewerSetFormat( view, PETSC_VIEWER_ASCII_SYMMODU );
    VecView( abs, view );

    VecDuplicate( *wf, &prob );
    VecAssemblyBegin( prob );
    VecAssemblyEnd( prob );

    std::vector<PetscReal> efvec;
    std::vector<std::ofstream> dipole;
    // dipol( 3 * cntx->dipole->decompositions().size() + 1 );
    if ( cntx->hparams->rank() == 0 ) {
        for ( auto a = 0; a < cntx->dipole->dipole_filename().size();
              ++a ) {
            try
            {
                dipole.emplace_back(
                    std::ofstream( cntx->dipole->dipole_filename()[a],
                                   std::ios::binary ) );
            }
            catch ( ... )
            {
                std::cerr << "couldn't open the dipole file, oops... it "
                             "won't be "
                             "written to disk" << std::endl;
            }
        }
    }
    // std::vector< std::vector<PetscReal> > dipole( 3 *
    // cntx->dipole->decompositions().size() + 1 );
    std::vector<PetscReal> time;

    PetscReal norm_lost = 0;

    while ( t < maxtime ) {
        MatCopy( *( cntx->D ), *A, SAME_NONZERO_PATTERN ); // A = D

        // ef-forward:
        PetscScalar eff = cntx->laser->efield( t + cntx->laser->dt() );
        // This has different 't's on both sides:
        MatScale( *A, ( ef + eff ) / 2. ); // A = ef(t) * D
        MatDiagonalSet(
            *A, *( cntx->H ), INSERT_VALUES ); // A = ef(t) * D + H_0
        MatScale( *A, cn_factor ); // A = - i * .5 * dt [ ef(t) * D + H_0 ]
        MatShift( *A,
                  std::complex<double>( 1, 0 ) ); // A = - i * .5 * dt [
                                                  // ef(t) * D + H_0 ] + 1
        MatMult( *A, *wf, tmp );                  // A u_n = tmp
        MatScale( *A, std::complex<double>( -1, 0 ) ); // A = i * .5 dt
        // [ef(t+dt) * D + H_0 ]
        // - 1
        MatShift(
            *A, std::complex<double>( 2, 0 ) ); // A = i * .5 dt [ef(t+dt)
                                                // * D + H_0 ] + 1


        KSPSetOperators( ksp, *A, *A, SAME_NONZERO_PATTERN ); // Solve[ A x
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
        t += cntx->laser->dt();
        ef = cntx->laser->efield( t );

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
            PetscViewerBinaryOpen( cntx->hparams->comm(),
                                   wf_name.str().c_str(),
                                   FILE_MODE_WRITE,
                                   &view );
            zero++;
            // analytically propagate for said time:
            std::cout << "analytically propagating the pulse: ";
            PetscReal dt = cntx->laser->next_pulse_start( t ) -
                           t; // the difference between then and now.
            std::cout << dt << std::endl;
            math::FieldFreePropagate(
                ( cntx->H ), wf, dt ); // propagate forward in time
            t = cntx->laser->next_pulse_start( t );
            ef = cntx->laser->efield( t );
            // write out the next zero:
            if ( cntx->hparams->rank() == 0 )
                std::cout << "time: " << t << " step: " << step
                          << " efield: " << ef << " norm-1: " << norm - 1
                          << " *" << zero << std::endl;
            wf_name.str( "" );
            wf_name << "./wf_" << zero << ".dat";
            PetscViewerBinaryOpen( cntx->hparams->comm(),
                                   wf_name.str().c_str(),
                                   FILE_MODE_WRITE,
                                   &view );
            // VecView( *wf, view );
            zero++;
            std::cout << output::reset;
            continue;
        }
        // else do nothing...

        step++;
        // at the zero of the field: write out the vector:
        if ( ( ef.real() <= 0. &&
               cntx->laser->efield( t - cntx->laser->dt() ).real() > 0 ) ||
             ( ef.real() >= 0 &&
               cntx->laser->efield( t - cntx->laser->dt() ).real() <
                   0 ) ) {
            if ( cntx->hparams->rank() == 0 )
                std::cout << "time: " << t << " step: " << step
                          << " efield: " << ef << " norm-1: " << norm - 1
                          << " *" << std::endl;
            std::ostringstream wf_name;
            wf_name << "./wf_" << zero << ".dat";
            PetscViewerBinaryOpen( cntx->hparams->comm(),
                                   wf_name.str().c_str(),
                                   FILE_MODE_WRITE,
                                   &view );
            VecView( *wf, view );
            zero++;
        }

        if ( true ) //(!(step%10))
        {
            if ( cntx->hparams->rank() == 0 ) {
                time.push_back( t );
                efvec.push_back( ef.real() );
            }
            cntx->dipole->find_dipole_moment_decompositions(
                *( cntx->D ), *wf, dipole, cntx->hparams->prototype() );
        }
        if ( !( step % 100 ) || step < 100 ) {
            VecCopy( *wf, prob );
            VecAbs( prob );
            VecNorm( prob, NORM_2, &norm );
            if ( cntx->hparams->rank() == 0 )
                std::cout << "time: " << t << " step: " << step
                          << " efield: " << ef << " norm-1: " << norm - 1
                          << " norm lost to absorbers: " << norm_lost
                          << std::endl;
        }
    }
    file_name = std::string( "./wf_final.dat" );
    PetscViewerBinaryOpen(
        cntx->hparams->comm(), file_name.c_str(), FILE_MODE_WRITE, &view );
    VecView( *wf, view );

    // we don't want to crash the nodes from lack of memory, so we only do
    // this on the head node
    if ( cntx->hparams->rank() == 0 ) {
        try
        {
            common::export_vector_binary( "time.dat", time );
        }
        catch ( ... )
        {
            std::cerr
                << "couldn't open the time file, oops... it won't be "
                   "written to disk" << std::endl;
        }
        try
        {
            common::export_vector_binary( cntx->laser->laser_filename(),
                                          efvec );
        }
        catch ( ... )
        {
            std::cerr
                << "couldn't open the efvec file, oops... it won't be "
                   "written to disk" << std::endl;
        }

        // close the dipole files
        for ( auto &a : dipole ) {
            a.close();
        }
    }

    KSPDestroy( &ksp );
    std::cerr << "leaving cranknicholson" << std::endl;
    return 0;
}
}
