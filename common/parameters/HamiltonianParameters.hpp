#pragma once

//ebss:
#include<common/common.hpp>
#include<common/parameters.hpp>

//petsc:
#include<petsc.h>

template < typename write_type_ >
class HamiltonianParameters: public Parameters
{
public:

    typedef struct {
        PetscInt nmax, lmax, mmax;
        char basis_folder[PETSC_MAX_PATH_LEN];
        char dipole_matrix_filename[PETSC_MAX_PATH_LEN];
        char energy_eigenvalues_filename[PETSC_MAX_PATH_LEN];
        char prototype_filename[PETSC_MAX_PATH_LEN];
        char basis_function_bag[PETSC_MAX_PATH_LEN];
    } hamiltonian;

    HamiltonianParameters(MPI_Comm comm): Parameters(comm)
    {
        this->params = new hamiltonian;
        p = (void**) params;
        PetscBagCreate(this->comm, sizeof(hamiltonian), &this->bag);
        PetscBagGetData(this->bag, this->p);
        this->register_params();
        PetscBagSetFromOptions(this->bag);
    }
    HamiltonianParameters(MPI_Comm comm, std::string filename): Parameters(comm)
    {
        this->params = new hamiltonian;
        p = (void**) params;
        PetscBagCreate(this->comm, sizeof(hamiltonian), &this->bag);
        PetscBagGetData(this->bag, this->p);
        this->register_params();
        PetscBagSetFromOptions(this->bag);
    }

    PetscErrorCode save_parameters();
    PetscErrorCode save_parameters(std::string filename);
    PetscErrorCode init_from_file();
    PetscErrorCode init_from_file(std::string filename);

    ~BasisParameters()
    {
        this->save_parameters();
        delete params;
    }

    PetscInt nmax() const;
    PetscInt lmax() const;
    PetscInt mmax() const;

    std::string dipole_matrix_filename() const;
    std::string energy_eigenvalues_filename() const;
    std::string prototype_filename() const;
    BasisParameters<write_type_> basis_parameters();


private:
    BasisParameters<write_type_> basis;
    std::vector<BasisID> *prototype;

}



