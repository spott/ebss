#pragma once

#include<array>

typedef struct {
    Mat *D;
    Vec *H;
    HamiltonianParameters<PetscReal> *hparams;
    PulsetrainParameters *laser;
    AbsorberParameters *absorber;
    DipoleParameters *dipole;
} context;

namespace cranknicholson
{

PetscErrorCode
solve(Vec *wf, context* cntx, Mat *A)
{
    PetscViewer view;
    KSP     ksp;
    PC      pc;
    KSPCreate(cntx->hparams->comm(), &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCJACOBI);
    KSPSetTolerances(ksp,1e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(ksp);

    PetscScalar     cn_factor = std::complex<double>(0, -0.5 * cntx->laser->dt());
    PetscReal       t = 0;
    PetscScalar     ef = cntx->laser->efield(t);
    PetscReal       maxtime = cntx->laser->max_time();
    PetscInt        step = 0;
    Vec             tmp;
    Vec             prob;
    PetscReal       norm;
    PetscInt        zero = 0;

    if (cntx->hparams->rank() == 0)
        std::cout << "maxtime: " << maxtime << " steps: " << maxtime/cntx->laser->dt() << " pulse length: " << cntx->laser->pulse_length() << std::endl;

    VecDuplicate(*wf, &tmp);
    VecAssemblyBegin(tmp);
    VecAssemblyEnd(tmp);

    Vec abs;
    VecDuplicate(*wf, &abs);
    cntx->absorber->absorb(&abs, cntx->hparams);
    //if (cntx->absorber->type() == "cx_rot" || cntx->absorber->type() == "cx_scale")
        //VecPointwiseMult(*(cntx->H), *(cntx->H), abs);

    std::string file_name = std::string("./absorber.dat");
    PetscViewerASCIIOpen(cntx->hparams->comm(),file_name.c_str(),&view);
    PetscViewerSetFormat(view, PETSC_VIEWER_ASCII_SYMMODU);
    VecView(abs,view);

    VecDuplicate(*wf, &prob);
    VecAssemblyBegin(prob);
    VecAssemblyEnd(prob);

    std::vector< PetscReal > efvec;
    std::vector< PetscReal > dipole;
    std::vector< PetscReal > time;
    
    PetscReal norm_lost = 0;

    while (t < maxtime)
    {

        MatCopy(*( cntx->D ), *A , SAME_NONZERO_PATTERN);   // A = D

		//ef-forward:
		PetscScalar eff = cntx->laser->efield( t+ cntx->laser->dt() );
        //This has different 't's on both sides:
        MatScale(*A, (ef + eff)/2.  );                                   // A = ef(t) * D
        MatDiagonalSet(*A, *(cntx->H), INSERT_VALUES);      // A = ef(t) * D + H_0
        MatScale(*A, cn_factor);                            // A = - i * .5 * dt [ ef(t) * D + H_0 ]
        MatShift(*A, std::complex<double>(1,0));            // A = - i * .5 * dt [ ef(t) * D + H_0 ] + 1
        MatMult(*A, *wf, tmp);                              // A u_n = tmp
       	//MatAXPY(*A, cn_factor * (cntx->laser->efield(t+cntx->laser->dt()) - ef), *( cntx->D ), SAME_NONZERO_PATTERN);
        MatScale(*A, std::complex<double>(-1,0));           // A = i * .5 dt [ef(t+dt) * D + H_0 ] - 1
        MatShift(*A, std::complex<double>(2,0));            // A = i * .5 dt [ef(t+dt) * D + H_0 ] + 1

        KSPSetOperators(ksp, *A, *A, SAME_NONZERO_PATTERN); // Solve[ A x = tmp ] for x

        KSPSolve(ksp, tmp, *wf);

        if (cntx->absorber->type() == "cosine")
        {
            PetscReal n = 0;
            PetscReal n2 = 0;
            VecNorm(*wf,NORM_2,&n);
            VecPointwiseMult(*wf, *wf, abs);
            VecNorm(*wf,NORM_2,&n2);
            norm_lost += n2-n;
	        //if (cntx->hparams->rank() == 0) std::cout << "lost " << n2-n << " norm, total = " << norm_lost << std::endl;
        }

        //look at the next point
        t += cntx->laser->dt();
        ef = cntx->laser->efield(t);
        
        //check if we are "in" a pulse now:
        if(!cntx->laser->in_pulse(t))
        {
            //first write out the zero:
            std::cout << output::red ;
            if (cntx->hparams->rank() == 0)
                std::cout << "time: " << t << " step: " << step << " efield: " << ef << " norm-1: " << norm-1 << " * " << zero << std::endl;
            std::ostringstream wf_name;
            wf_name << "./wf_" << zero << ".dat";
            PetscViewerBinaryOpen(cntx->hparams->comm(),wf_name.str().c_str(),FILE_MODE_WRITE,&view);
            VecView(*wf, PETSC_VIEWER_STDOUT_WORLD);
            //PetscViewerSetFormat(view, PETSC_VIEWER_ASCII_MATHEMATICA);
            VecView(*wf, view);
            zero++;
            //analytically propagate for said time:
            std::cout << "analytically propagating the pulse: ";
            PetscReal dt = cntx->laser->next_pulse_start(t) - t; //the difference between then and now.
            std::cout << dt <<  std::endl;
            math::FieldFreePropagate((cntx->H), wf, dt); //propagate forward in time
            t = cntx->laser->next_pulse_start(t);
            ef = cntx->laser->efield(t);
            //write out the next zero:
            if (cntx->hparams->rank() == 0)
                std::cout << "time: " << t << " step: " << step << " efield: " << ef << " norm-1: " << norm-1 << " *" << zero << std::endl;
            wf_name.str("");
            wf_name << "./wf_" << zero << ".dat";
            PetscViewerBinaryOpen(cntx->hparams->comm(),wf_name.str().c_str(),FILE_MODE_WRITE,&view);
            //PetscViewerSetFormat(view, PETSC_VIEWER_ASCII_SYMMODU);
            VecView(*wf, view);
            zero++;
            std::cout << output::reset;
            continue;
        }
        //else do nothing...

        step++;
        //at the zero of the field: write out the vector:
        if (( ef.real() <= 0. && cntx->laser->efield(t - cntx->laser->dt()).real() > 0 ) || ( ef.real() >= 0 && cntx->laser->efield(t - cntx->laser->dt()).real() < 0 ))
        {
            if (cntx->hparams->rank() == 0)
                std::cout << "time: " << t << " step: " << step << " efield: " << ef << " norm-1: " << norm-1 << " *" << std::endl;
            std::ostringstream wf_name;
            wf_name << "./wf_" << zero << ".dat";
            PetscViewerBinaryOpen(cntx->hparams->comm(),wf_name.str().c_str(),FILE_MODE_WRITE,&view);
            //PetscViewerASCIIOpen(cntx->hparams->comm(),wf_name.str().c_str(),&view);
            //PetscViewerSetFormat(view, PETSC_VIEWER_ASCII_SYMMODU);
            VecView(*wf, view);
            zero++;
        }
        if (!(step%10))
        {
            time.push_back( t );
            efvec.push_back( ef.real() );
            dipole.push_back( cntx->dipole->find_dipole_moment(*(cntx->D), *wf) );
        }
        if (!(step%100))
        {
            VecCopy(*wf, prob);
            VecAbs(prob);
            VecNorm(prob,NORM_2,&norm);
            VecPointwiseMult(prob, prob, prob);
            VecShift(prob, 1e-20);
            VecLog(prob);
            //VecView(prob, PETSC_VIEWER_DRAW_WORLD);
            if (cntx->hparams->rank() == 0)
                std::cout << "time: " << t << " step: " << step << " efield: " << ef << " norm-1: " << norm-1 << " norm lost to absorbers: " <<  norm_lost << std::endl;
            //if (norm-1 > 10e-5 && cntx->hparams->rank() == 0)
                //std::cerr << "time: " << t << " step: " << step << " efield: " << ef << " norm-1: " << norm-1 << std::endl;
        }
    }
    file_name = std::string("./wf_final.dat");
    PetscViewerBinaryOpen(cntx->hparams->comm(),file_name.c_str(),FILE_MODE_WRITE,&view);
    VecView(*wf,view);

    if (cntx->hparams->rank() == 0) 
    {
        try {
            common::export_vector_binary( "time.dat" , time );
        } catch (...) {
            std::cerr << "couldn't open the time file, oops... it won't be written to disk" << std::endl;
        }
        try {
            common::export_vector_binary( cntx->laser->laser_filename() , efvec);
        } catch (...) {
            std::cerr << "couldn't open the efvec file, oops... it won't be written to disk" << std::endl;
        }

        try {
            common::export_vector_binary( cntx->dipole->dipole_filename() , dipole);
        } catch (...) {
            std::cerr << "couldn't open the dipole file, oops... it won't be written to disk" << std::endl;
        }
    }


    //After the propagation through diffeq, we need to do the propagation after the fact to find the dipole 
    //moment over time.
    

    std::vector< PetscReal > after_dipole;
    std::vector< PetscReal > after_time;
    while ( (t - maxtime) <= cntx->dipole->t_after() )
    {
        VecCopy(*wf, tmp);
        math::FieldFreePropagate(cntx->H, &tmp, t);
        after_time.push_back( t );
        after_dipole.push_back(  cntx->dipole->find_dipole_moment(*(cntx->D), tmp) );
        t += cntx->dipole->dt();
    }

    if (cntx->hparams->rank() == 0) 
    {
        try {
            common::export_vector_binary( cntx->dipole->after_filename() , after_dipole);
        } catch (...) {
            std::cerr << "couldn't open the after_dipole file, oops... it won't be written to disk... NOOO!!" << std::endl;
        }
        
        try {
            common::export_vector_binary( "after_time.dat" , after_time );
        } catch (...) {
            std::cerr << "couldn't open the after_time file, oops... it won't be written to disk... NOOO!!" << std::endl;
        }
    }

    //fourier transform and output the dipole timeseries:
    if (cntx->hparams->rank() == 0)
    {
        try {
            //this modifies the old vectors!  warning!
            common::export_vector_binary( 
                    cntx->dipole->dipole_filename()
                        .substr(0, cntx->dipole->dipole_filename().size() - 3)
                            .append("_fourier.dat") , 
                    math::fourier( math::window::hann_window(dipole ) ) );
            common::export_vector_binary( 
                    cntx->dipole->after_filename()
                        .substr(0, cntx->dipole->after_filename().size() - 3)
                            .append("_fourier.dat") , 
                    math::fourier( math::window::hann_window(after_dipole ) ) );

        } catch (...) {
            std::cerr << "couldn't write out the dipole_fourier files.. :( " << std::endl;
        }
    }

    KSPDestroy(&ksp);
    std::cerr << "leaving cranknicholson" << std::endl;
    return 0;
}


}
