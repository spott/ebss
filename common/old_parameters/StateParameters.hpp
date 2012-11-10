#pragma once

//ebss:
#include<common/common.hpp>
#include<common/math.hpp>
#include<common/parameters/Parameters.hpp>

//stl:
#include<sstream>
#include<string>

//petsc:
#include<petsc.h>

class StateParameters: public Parameters
{
public:
    typedef struct {
        PetscBool bound_states;
        PetscInt n[30];
        PetscInt l[30];
    } state;

    StateParameters(MPI_Comm comm): Parameters(comm)
    {
        PetscBool flg = PETSC_FALSE;
        char bagname[PETSC_MAX_PATH_LEN];
        PetscOptionsGetString(PETSC_NULL, "-state_bag_filename", bagname, PETSC_MAX_PATH_LEN, &flg);

        this->params = new state;

        if (flg)
        {
            this->bag_filename = std::string(bagname);
            this->init_from_file();
        }
        else
        {
            PetscBagCreate(this->comm_, sizeof(state), &this->bag);
            PetscBagGetData(this->bag, (void**)&params);
            this->register_params();
            PetscBagSetFromOptions(this->bag);
            //see if the folder exists:
            this->bag_filename = std::string("./State.bag");
        }
    }

    StateParameters(MPI_Comm comm, std::string filename): Parameters(comm)
    {
    }

    PetscInt n_size() const;
    PetscInt l_size() const;
    PetscInt m_size() const;
    PetscReal cos_factor() const;

    PetscErrorCode init_from_file();
    PetscErrorCode save_parameters();

private:
    state* params;
    PetscErrorCode register_params();
};


PetscErrorCode StateParameters::init_from_file()
{
    PetscErrorCode e = PetscBagCreate(this->comm_, sizeof(state), &this->bag);
    e = PetscBagGetData(this->bag, (void**)&params);
    this->register_params();

    e = Parameters::init_from_file(this->bag_filename);
    return e;
}

PetscInt StateParameters::n_size() const { return params->n_size; }
PetscInt StateParameters::l_size() const { return params->l_size; }
PetscInt StateParameters::m_size() const { return params->m_size; }
PetscReal StateParameters::cos_factor() const { return params->cos_factor; }

PetscErrorCode StateParameters::save_parameters()
{
    return Parameters::save_parameters(this->bag_filename);
}

PetscErrorCode 
StateParameters::register_params()
{
    PetscErrorCode ierr;
    ierr = PetscBagSetOptionsPrefix(this->bag, "state_" );
    ierr = PetscBagSetName(this->bag,
            "StateParams",
            "Parameters for modifying the dipole transitions");
    ierr = PetscBagRegisterInt(this->bag, &params->n_size,
            40, "n_size", "size of the absorber in n");
    ierr = PetscBagRegisterInt(this->bag, &params->l_size,
            10, "l_size", "size of the absorber in l");
    ierr = PetscBagRegisterInt(this->bag, &params->m_size,
            0, "m_size", "size of the absorber in m");
    ierr = PetscBagRegisterReal(this->bag, &params->cos_factor,
            0.125, "cos_factor", "the exponential on the cosine");
    return ierr;
}
