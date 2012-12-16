#pragma once

typedef struct {
    Mat *D;
    Vec *H;
    HamiltonianParameters<PetscReal> *hparams;
    PulsetrainParameters *laser;
    AbsorberParameters *absorber;
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

    PetscScalar     cn_factor = std::complex<double>(0, -0.5 * cntx->laser->dt());
    PetscReal       t = 0;
    PetscScalar     ef = cntx->laser->efield(t); // + cntx->laser->efield(t+cntx->laser->cycles());  //average between the points
    //ef /= 2;
    PetscReal       maxtime = cntx->laser->max_time();
    PetscInt        step = 0;
    Vec             tmp;
    Vec             prob;
    PetscReal       norm;
    PetscInt        zero = 0;

    if (cntx->hparams->rank() == 0)
        std::cout << "maxtime: " << maxtime << " steps: " << maxtime/cntx->laser->dt() << std::endl;

    VecDuplicate(*wf, &tmp);
    VecAssemblyBegin(tmp);
    VecAssemblyEnd(tmp);

    Vec abs;
    VecDuplicate(*wf, &abs);
    cntx->absorber->absorb(&abs, cntx->hparams);
    if (cntx->absorber->type() == "cx_rot" || cntx->absorber->type() == "cx_scale")
        VecPointwiseMult(*(cntx->H), *(cntx->H), abs);

    std::string file_name = std::string("./absorber.dat");
    PetscViewerASCIIOpen(cntx->hparams->comm(),file_name.c_str(),&view);
    PetscViewerSetFormat(view, PETSC_VIEWER_ASCII_SYMMODU);
    VecView(abs,view);

    VecDuplicate(*wf, &prob);
    VecAssemblyBegin(prob);
    VecAssemblyEnd(prob);

    std::vector<PetscReal> efvec;

    while (t < maxtime)
    {
        efvec.push_back(ef.real());
        MatCopy(*( cntx->D ), *A , SAME_NONZERO_PATTERN);   // A = D

        //This has the same "t" at both sides of the equation... should be different...
        //MatScale(*A, ef);                                   // A = ef(t) * D
        //MatDiagonalSet(*A, *(cntx->H), INSERT_VALUES);      // A = ef(t) * D + H_0
        //MatScale(*A, cn_factor);                            // A = .5 * dt [ ef(t) * D + H_0 ]
        //MatShift(*A, std::complex<double>(1,0));            // A = .5 * dt [ ef(t) * D + H_0 ] + 1
        //MatMult(*A, *wf, tmp);                              // A u_n = tmp
        //MatScale(*A, std::complex<double>(-1,0));           // A = - .5 dt [ef(t+dt) * D + H_0 ] - 1
        //MatShift(*A, std::complex<double>(2,0));            // A = - .5 dt [ef(t+dt) * D + H_0 ] + 1
        
        //This has different 't's on both sides:
        MatScale(*A, ef);                                   // A = ef(t) * D
        MatDiagonalSet(*A, *(cntx->H), INSERT_VALUES);      // A = ef(t) * D + H_0
        MatScale(*A, cn_factor);                            // A = .5 * dt [ ef(t) * D + H_0 ]
        MatShift(*A, std::complex<double>(1,0));            // A = .5 * dt [ ef(t) * D + H_0 ] + 1
        MatMult(*A, *wf, tmp);                              // A u_n = tmp
        MatAXPY(*A, cn_factor * (cntx->laser->efield(t+cntx->laser->dt()) - ef), *( cntx->D ), SAME_NONZERO_PATTERN);
        MatScale(*A, std::complex<double>(-1,0));           // A = - .5 dt [ef(t+dt) * D + H_0 ] - 1
        MatShift(*A, std::complex<double>(2,0));            // A = - .5 dt [ef(t+dt) * D + H_0 ] + 1

        KSPSetOperators(ksp, *A, *A, SAME_NONZERO_PATTERN); // Solve[ A x = tmp ] for x
        KSPSetFromOptions(ksp);

        KSPSolve(ksp, tmp, *wf);

        if (cntx->absorber->type() == "cosine")
            VecPointwiseMult(*wf, *wf, abs);

        t += cntx->laser->dt();
        ef = cntx->laser->efield(t); //  + cntx->laser->efield(t+cntx->laser->cycles());  //average between the points
        //ef /= 2;
        step++;
        //at the zero of the field: write out the vector:
        if (( ef.real() <= 0. && cntx->laser->efield(t - cntx->laser->dt()).real() > 0 ) || ( ef.real() >= 0 && cntx->laser->efield(t - cntx->laser->dt()).real() < 0 ))
        {
            if (cntx->hparams->rank() == 0)
                std::cout << "time: " << t << " step: " << step << " efield: " << ef << " norm-1: " << norm-1 << " *" << std::endl;
            std::ostringstream wf_name;
            wf_name << "./wf_" << zero << ".dat";
            PetscViewerASCIIOpen(cntx->hparams->comm(),wf_name.str().c_str(),&view);
            PetscViewerSetFormat(view, PETSC_VIEWER_ASCII_SYMMODU);
            VecView(*wf, view);
            zero++;
        }

        if (!(step%50))
        {
            VecCopy(*wf, prob);
            VecAbs(prob);
            VecNorm(prob,NORM_2,&norm);
            VecPointwiseMult(prob, prob, prob);
            VecShift(prob, 1e-20);
            VecLog(prob);
            VecView(prob, PETSC_VIEWER_DRAW_WORLD);
            if (cntx->hparams->rank() == 0)
                std::cout << "time: " << t << " step: " << step << " efield: " << ef << " norm-1: " << norm-1 << std::endl;
            //if (norm-1 > 10e-5 && cntx->hparams->rank() == 0)
                //std::cerr << "time: " << t << " step: " << step << " efield: " << ef << " norm-1: " << norm-1 << std::endl;
        }
    }
    file_name = std::string("./wf_final.dat");
    PetscViewerASCIIOpen(cntx->hparams->comm(),file_name.c_str(),&view);
    PetscViewerSetFormat(view, PETSC_VIEWER_ASCII_SYMMODU);
    VecView(*wf,view);

    if (cntx->hparams->rank() == 0) common::export_vector_ascii( cntx->laser->laser_filename() , efvec);
    KSPDestroy(&ksp);
    std::cerr << "leaving cranknicholson" << std::endl;
    return 0;
}

}
