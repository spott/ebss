#pragma once

#include<petsc.h>


class Parameters
{
public:
    Parameters(MPI_Comm comm): comm_(comm) {};
    Parameters(MPI_Comm comm, void** s): comm_(comm), p(s)  {};

    PetscErrorCode init_from_file(const char* filename);
    PetscErrorCode save_parameters(const char* filename);
    MPI_Comm comm() const;
    PetscErrorCode print_parameters();

    virtual ~Parameters() {};

protected:
    virtual PetscErrorCode register_params() { return 0; };
    MPI_Comm        comm_;
    PetscBag        bag;
    void**          p;
};

PetscErrorCode  Parameters::init_from_file(const char* filename)
{
    PetscViewer viewer;
    PetscViewerBinaryOpen(this->comm_,filename,FILE_MODE_READ,&viewer);
    PetscBagLoad(viewer, &this->bag);
    PetscBagGetData(this->bag, (void **)(p));
    PetscBagSetFromOptions(this->bag);
    PetscViewerDestroy(&viewer);
    return 0;
};

PetscErrorCode  Parameters::save_parameters(const char* filename)
{
    PetscErrorCode ierr;
    PetscViewer viewer;
    ierr = PetscViewerBinaryOpen(this->comm_,filename,FILE_MODE_WRITE,&viewer);
    ierr = PetscBagView(this->bag, viewer);
    ierr = PetscViewerDestroy(&viewer);
    return ierr;
};

MPI_Comm Parameters::comm() const 
{ 
    return comm_; 
};

PetscErrorCode Parameters::print_parameters()
{
    PetscErrorCode ierr;
    ierr = PetscBagView(this->bag,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    return (0);
};


class BasisParameters: public Parameters
{
public:
    typedef struct {
        PetscReal rmax, rmin;
        PetscInt  lmax, nmax, points;
        char basis_folder[PETSC_MAX_PATH_LEN];
    } basis;

    BasisParameters(MPI_Comm comm): Parameters(comm)
    {
        this->params = new basis;
        p = (void**) params;
        PetscBagCreate(this->comm_, sizeof(basis), &this->bag);
        PetscBagGetData(this->bag, this->p);
        this->register_params();
        PetscBagSetFromOptions(this->bag);
    };

private:
    //std::vector<PetscScalar> grid;
    PetscErrorCode register_params();
    basis* params;

};
    

PetscErrorCode BasisParameters::register_params()
{
    PetscErrorCode ierr;
    ierr = PetscBagSetName(this->bag,
            "BasisParams",
            "Parameters for finding the basis state");
    ierr = PetscBagRegisterInt(this->bag, &params->nmax,
            500, "nmax", "Max n value");
    ierr = PetscBagRegisterInt(this->bag, &params->lmax,
            50, "lmax", "Max l value");
    ierr = PetscBagRegisterReal(this->bag, &params->rmax,
            1000., "rmax", "Maximum r for grid");
    ierr = PetscBagRegisterReal(this->bag, &params->rmin,
            .1, "rmin", "Minimum r for grid");
    ierr = PetscBagRegisterInt(this->bag, &params->points,
            10000, "points", "Number of points");
    ierr = PetscBagRegisterString(this->bag, &params->basis_folder,
            PETSC_MAX_PATH_LEN, "./basis/",
            "evectors_folder",
            "Where the evectors and evalues should be stored");
    return ierr;
}
