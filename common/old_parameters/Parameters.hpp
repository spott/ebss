#pragma once

#include<iostream>
#include<string>
#include<algorithm>
#include<petsc.h>

struct BasisID{
    PetscInt n;
    PetscInt l;
    PetscInt m;
    PetscScalar e;
    bool operator<(const BasisID b)
    {
        if (this->l < b.l)
            return true;
        else if (this->l == b.l && this->n < b.n)
            return true;
        else if (this->l == b.l && this->n == b.n && this->m < b.m)
            return true;
        else
            return false;
    }
};
bool operator!=(const BasisID &a, const BasisID &b)
{
    if (a.l != b.l || a.n != b.n || a.m != b.m || a.e != b.e)
        return true;
    else
        return false;
}
std::ostream& operator<<(std::ostream &out, const BasisID &b)     //output
{
    out << b.n << ", " << b.l << ", " << b.m << ", " << b.e;
    return out;
}

class Parameters
{
public:
    Parameters(MPI_Comm comm): comm_(comm) {};

    virtual PetscErrorCode init_from_file(std::string filename);
    virtual PetscErrorCode save_parameters(std::string filename);
    MPI_Comm comm() const;
    int rank() const;
    int size() const;
    virtual PetscErrorCode print_parameters();

    virtual ~Parameters() {};

protected:
    virtual PetscErrorCode register_params() { return 0; };
    std::string     bag_filename;
    MPI_Comm        comm_;
    PetscBag        bag;
    //void**          p;
};

PetscErrorCode  Parameters::init_from_file(std::string filename)
{
    PetscViewer viewer;
    PetscViewerBinaryOpen(this->comm_,filename.c_str(),FILE_MODE_READ,&viewer);
    PetscBagLoad(viewer, this->bag);
    //this->register_params();
    //PetscBagSetFromOptions(this->bag);
    PetscViewerDestroy(&viewer);
    return 0;
};

PetscErrorCode  Parameters::save_parameters(std::string filename)
{
    PetscErrorCode ierr;
    PetscViewer viewer;
    ierr = PetscViewerBinaryOpen(this->comm_,filename.c_str(),FILE_MODE_WRITE,&viewer);
    ierr = PetscBagView(this->bag, viewer);
    ierr = PetscViewerDestroy(&viewer);
    return ierr;
};

int Parameters::rank() const
{
    int i = -1;
    MPI_Comm_rank(this->comm_, &i);

    if (i < 0)
        throw (std::exception());
    return i;
}

int Parameters::size() const
{
    int i = -1;
    MPI_Comm_size(this->comm_, &i);

    if (i < 0)
        throw (std::exception());
    return i;
}

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

