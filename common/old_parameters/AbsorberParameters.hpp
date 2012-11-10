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

class AbsorberParameters: public Parameters
{
public:
    typedef struct {
        PetscInt n_size;
        PetscInt l_size;
        PetscInt m_size;
        PetscReal cos_factor;
    } absorber;

    AbsorberParameters(MPI_Comm comm): Parameters(comm)
    {
        PetscBool flg = PETSC_FALSE;
        char bagname[PETSC_MAX_PATH_LEN];
        PetscOptionsGetString(PETSC_NULL, "-absorber_bag_filename", bagname, PETSC_MAX_PATH_LEN, &flg);

        this->params = new absorber;

        if (flg)
        {
            this->bag_filename = std::string(bagname);
            this->init_from_file();
        }
        else
        {
            PetscBagCreate(this->comm_, sizeof(absorber), &this->bag);
            PetscBagGetData(this->bag, (void**)&params);
            this->register_params();
            PetscBagSetFromOptions(this->bag);
            //see if the folder exists:
            this->bag_filename = std::string("./Absorber.bag");
        }
    }

    AbsorberParameters(MPI_Comm comm, std::string filename): Parameters(comm)
    {
    }

    PetscInt n_size() const;
    PetscInt l_size() const;
    PetscInt m_size() const;
    PetscReal cos_factor() const;

    PetscErrorCode init_from_file();
    PetscErrorCode save_parameters();

private:
    absorber* params;
    PetscErrorCode register_params();
};


PetscErrorCode AbsorberParameters::init_from_file()
{
    PetscErrorCode e = PetscBagCreate(this->comm_, sizeof(absorber), &this->bag);
    e = PetscBagGetData(this->bag, (void**)&params);
    this->register_params();

    e = Parameters::init_from_file(this->bag_filename);
    return e;
}

PetscInt AbsorberParameters::n_size() const { return params->n_size; }
PetscInt AbsorberParameters::l_size() const { return params->l_size; }
PetscInt AbsorberParameters::m_size() const { return params->m_size; }
PetscReal AbsorberParameters::cos_factor() const { return params->cos_factor; }

PetscErrorCode AbsorberParameters::save_parameters()
{
    return Parameters::save_parameters(this->bag_filename);
}

PetscErrorCode 
AbsorberParameters::register_params()
{
    PetscErrorCode ierr;
    ierr = PetscBagSetOptionsPrefix(this->bag, "absorber_" );
    ierr = PetscBagSetName(this->bag,
            "AbsorberParams",
            "Parameters for the absorber");
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
