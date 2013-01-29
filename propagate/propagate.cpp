
//ebss:
#include<common/parameters/HamiltonianParameters.hpp>
#include<common/parameters/PulsetrainParameters.hpp>
#include<common/parameters/StateParameters.hpp>
#include<common/parameters/AbsorberParameters.hpp>
#include<common/parameters/DipoleParameters.hpp>
#include<common/common.hpp>
#include<common/output.hpp>
#include<propagate/cranknicholson.hpp>

//petsc:
#include<petsc.h>

//stl:
#include<sstream>
#include<string>
#include<iostream>


PetscErrorCode
Monitor(TS ts, PetscInt steps, PetscReal time, Vec x,void *ctx);

PetscErrorCode
HamiltonianJ(TS ts, PetscReal t, Vec u, Mat *A, Mat *B, MatStructure *flag, void *ctx);

int
main(int argc, const char ** argv)
{
    int ac = argc;
    char** av = new char*[argc];
    for (size_t i = 0; i < argc; i++)
    {
        av[i] = new char[strlen(argv[i])+1];
        std::copy(argv[i], argv[i] + strlen(argv[i])+1, av[i]);
    }
    PetscInitialize(&ac, &av, PETSC_NULL, PETSC_NULL);

    PetscViewer view;
    PetscBool flg = PETSC_FALSE;
    char bagname[PETSC_MAX_PATH_LEN];
    PetscOptionsGetString(PETSC_NULL, "-hamiltonian_config", bagname, PETSC_MAX_PATH_LEN, &flg);

    if (!flg)
    {
        std::cerr << "I need a hamiltonian to propagate. (-hamiltonian_config )" << std::endl;
        PetscFinalize();
        return 0;
    }

    HamiltonianParameters<PetscReal> *params = new HamiltonianParameters<PetscReal>(MPI_COMM_WORLD, std::string(bagname) );
    PulsetrainParameters *lparams = new PulsetrainParameters(argc, argv, MPI_COMM_WORLD);
    AbsorberParameters *aparams = new AbsorberParameters(argc, argv, MPI_COMM_WORLD);
    StateParameters *sparams = new StateParameters(argc, argv, MPI_COMM_WORLD);
    DipoleParameters *dparams = new DipoleParameters(argc, argv, MPI_COMM_WORLD);

    auto empty_states_index  = sparams->empty_states_index( params->prototype() );
    if (params->rank() == 0)
    {
        std::cout << params->print();
        std::cout << lparams->print();
        std::cout << aparams->print();
        std::cout << sparams->print();
        std::cout << dparams->print();

        common::export_vector_ascii(std::string("./prototype.csv"), params->prototype() );
        //params->save_parameters();
        lparams->save_parameters();
        aparams->save_parameters();
        sparams->save_parameters();
        dparams->save_parameters();
    }


    Mat D; //The Dipole Matrix
    Vec H; //The eigenvalues of the field free hamiltonian
    Vec wf; //The vec we are going to propagate
    Mat A; //The matrix we are using durring the solving
    TS ts; //The timestep context
    context* cntx = new context; //The context that we will pass around

    //fill the context:
    cntx->hparams = params;
    cntx->laser = lparams;
    cntx->absorber = aparams;
    cntx->dipole = dparams;
    cntx->H = &H;
    cntx->D = &D;

    D = params->read_dipole_matrix();
    MatAssemblyBegin(D, MAT_FINAL_ASSEMBLY);
    H = params->read_energy_eigenvalues();
    VecAssemblyBegin(H);

    MatAssemblyEnd(D,MAT_FINAL_ASSEMBLY);
    VecAssemblyEnd(H);

    //Do the state stuff... remove rows/columns:
    MatSetOption(D,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);
    MatZeroRowsColumns(D, empty_states_index.size(), empty_states_index.data(), 0.0, PETSC_NULL, PETSC_NULL);
    std::vector<PetscScalar> zeros(empty_states_index.size(), 0.0);
    if (zeros.size() != 0)
        VecSetValues(H, empty_states_index.size(), empty_states_index.data(), zeros.data(), INSERT_VALUES);


    VecAssemblyBegin(H);
    MatAssemblyBegin(D, MAT_FINAL_ASSEMBLY);
    VecAssemblyEnd(H);
    MatAssemblyEnd(D, MAT_FINAL_ASSEMBLY);
    std::cout << output::red << "H: " << std::endl;
    VecView(H, PETSC_VIEWER_STDOUT_WORLD);
    std::cout << output::reset <<std::endl;
    MatSetOption(D,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);

    //Setup the wavefunction:
    MatGetVecs(D, &wf, PETSC_NULL);
    sparams->initial_vector(&wf, params->prototype());
    std::cout << output::blue << "WF: " << std::endl;
    VecView(wf, PETSC_VIEWER_STDOUT_WORLD);
    std::cout << output::reset <<std::endl;


    //Copy the non-zero pattern from D to A (the matrix we use in our solver)
    MatDuplicate(D,MAT_SHARE_NONZERO_PATTERN,&A);
    MatSetFromOptions(A);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    //list where the eigenvalues are:
    Vec gg = common::eigen_balls(D);
    std::string file_name = std::string("./balls.dat");
    PetscViewerASCIIOpen(MPI_COMM_WORLD,file_name.c_str(),&view);
    PetscViewerSetFormat(view, PETSC_VIEWER_ASCII_SYMMODU);
    VecView(gg,view);

    //Propagate:

    cranknicholson::solve(&wf, cntx, &A);

    //create a viewer in the current directory:
    //file_name = std::string("./final_wf.dat");
    //PetscViewerASCIIOpen(MPI_COMM_WORLD,file_name.c_str(),&view);
    //PetscViewerSetFormat(view, PETSC_VIEWER_ASCII_SYMMODU);
    //VecView(wf,view);

    ////Create the TS context:
    //TSCreate(params->comm(), &ts);
    ////make it linear
    //TSSetProblemType(ts, TS_LINEAR);
    ////attach my monitor
    //TSMonitorSet(ts, Monitor, cntx, PETSC_NULL);
    ////make it a crank-nicholson solver:
    //TSSetType(ts, TSCN);

    ////the RHS is the jacobian because it linearly depends on wf
    //TSSetRHSFunction(ts, PETSC_NULL, TSComputeRHSFunctionLinear, cntx);
    ////The jacobian is where we do our solutions stuff
    //TSSetRHSJacobian(ts, A, A, HamiltonianJ, cntx);

    //TSSetInitialTimeStep(ts, 0.0, lparams->dt());
    //TSSetSolution(ts, wf);

    //VecSetValue(wf, 0, 1., INSERT_VALUES);
    //VecAssemblyBegin(wf);
    //VecAssemblyEnd(wf);

    ////TSSetDuration(ts, 5, 5 * .01);
    //TSSetDuration(ts, 
            //math::PI * lparams->cycles() / (lparams->frequency() * lparams->dt()),
            //math::PI * lparams->cycles() / lparams->frequency());
    //TSSetFromOptions(ts);

    //PetscReal time;
    //PetscInt steps;
    //TSSolve(ts, wf, &time);
    //TSGetTimeStepNumber(ts, &steps);


    delete params;
    PetscFinalize();
    return 0;
}

PetscErrorCode
Monitor(TS ts, PetscInt steps, PetscReal time, Vec x, void *ctx)
{
    context* cntx = (context*)ctx;
    PetscReal norm;

    VecNorm(x,NORM_2,&norm);
    if (cntx->hparams->rank() == 0)
        std::cout << "time: " << time << " step: " << steps << " norm-1: " << norm-1 << std::endl;

    VecView(x, PETSC_VIEWER_DRAW_WORLD);
    return 0;
}

PetscErrorCode
HamiltonianJ(TS ts, PetscReal t, Vec u, Mat *A, Mat *B, MatStructure *flag, void *ctx)
{
    Mat AA = *A;
    PetscErrorCode err;
    context* cntx = (context*)ctx;
    PetscScalar ef = cntx->laser->efield(t);
    if (cntx->hparams->rank() == 0)
        std::cout << "ef: " << ef << std::endl;

    MatCopy(*(cntx->D), AA, SAME_NONZERO_PATTERN);
    MatScale(AA, ef);
    //MatZeroEntries(AA);
	err = MatDiagonalSet(AA,*(cntx->H),INSERT_VALUES);
	//err = MatAXPY(AA, ef, *(cntx->D), SAME_NONZERO_PATTERN);
    *flag = SAME_NONZERO_PATTERN;

    //A = B;
    if (cntx->hparams->rank() == 0)
        std::cout << "ef: " << ef << std::endl;

    return err;
}
