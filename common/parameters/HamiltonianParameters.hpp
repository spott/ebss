#pragma once

//ebss:
#include<common/common.hpp>
//#include<common/parameters/Parameters.hpp>
#include<common/parameters/BasisParameters.hpp>

//stl
#include<vector>
#include<sstream>
#include<string>

//petsc:
#include<petsc.h>

template < typename write_type_ >
class HamiltonianParameters: public Parameters
{
public:

    typedef struct {
        PetscInt nmax, lmax, mmax;
        char hamiltonian_folder[PETSC_MAX_PATH_LEN];
        char basis_bag[PETSC_MAX_PATH_LEN];
    } hamiltonian;

    HamiltonianParameters(MPI_Comm comm): Parameters(comm)
    {
        PetscBool flg = PETSC_FALSE;
        char bagname[PETSC_MAX_PATH_LEN];
        PetscOptionsGetString(PETSC_NULL, "-hamiltonian_bag_filename", bagname, PETSC_MAX_PATH_LEN, &flg);

        this->params = new hamiltonian;
        p = (void**) params;
        if (flg)
        {
            this->bag_filename = std::string(bagname);
            this->init_from_file();
            this->basis_ = new BasisParameters<write_type_,write_type_>(comm, std::string(this->params->basis_bag));
        }
        else
        {
            PetscBagCreate(this->comm_, sizeof(hamiltonian), &this->bag);
            PetscBagGetData(this->bag, this->p);
            this->register_params();
            PetscBagSetFromOptions(this->bag);
            //see if the folder exists:
            boost::filesystem::create_directories(this->hamiltonian_folder());
            PetscOptionsGetString(PETSC_NULL, "-basis_bag_filename", bagname, PETSC_MAX_PATH_LEN, &flg);
            if (flg)
                this->basis_ = new BasisParameters<write_type_,write_type_>(comm);
            else
                this->basis_ = new BasisParameters<write_type_,write_type_>(comm, std::string(this->params->basis_bag));
        }
    }

    HamiltonianParameters(MPI_Comm comm, BasisParameters<write_type_>* basis): Parameters(comm), basis_(basis)
    {
        this->params = new hamiltonian;
        p = (void**) params;
        PetscBagCreate(this->comm_, sizeof(hamiltonian), &this->bag);
        PetscBagGetData(this->bag, this->p);
        this->register_params();
        PetscBagSetFromOptions(this->bag);
    }

    HamiltonianParameters(MPI_Comm comm, BasisParameters<write_type_>* basis, std::string filename): Parameters(comm), basis_(basis)
    {
        this->params = new hamiltonian;
        p = (void**) params;
        PetscBagCreate(this->comm, sizeof(hamiltonian), &this->bag);
        PetscBagGetData(this->bag, this->p);
        this->register_params();
        PetscBagSetFromOptions(this->bag);
    }

    PetscErrorCode save_parameters();
    PetscErrorCode init_from_file();

    ~HamiltonianParameters()
    {
        this->save_parameters();
        delete params;
    }

    PetscInt nmax() const;
    PetscInt lmax() const;
    PetscInt mmax() const;

    std::string hamiltonian_folder() const;
    std::string dipole_matrix_filename() const;
    std::string energy_eigenvalues_filename() const;
    std::string prototype_filename() const;
    BasisParameters<write_type_, write_type_> basis_parameters();

    std::vector<BasisID>* prototype() const;


private:
    BasisParameters<write_type_, write_type_> *basis_;
    std::vector<BasisID> *prototype_;
    hamiltonian* params;
    PetscErrorCode register_params();
};

template< typename write_type_ >
std::string HamiltonianParameters<write_type_>::hamiltonian_folder() const {return std::string(this->params->hamiltonian_folder);}
template< typename write_type_ >
std::string HamiltonianParameters<write_type_>::dipole_matrix_filename() const
{
    std::stringstream ss;
    std::string s;
    ss << params->hamiltonian_folder;
    ss << "dipole_matrix.dat";
    ss >> s;
    return s;
}


template< typename write_type_ >
std::string HamiltonianParameters<write_type_>::energy_eigenvalues_filename() const
{
    std::stringstream ss;
    std::string s;
    ss << params->hamiltonian_folder;
    ss << "energy_eigenvalues_vector.dat";
    ss >> s;
    return s;
}


template< typename write_type_ >
std::string HamiltonianParameters<write_type_>::prototype_filename() const
{
    std::stringstream ss;
    std::string s;
    ss << params->hamiltonian_folder;
    ss << "vector_prototype.dat";
    ss >> s;
    return s;
}

template< typename write_type_ >
BasisParameters<write_type_, write_type_> HamiltonianParameters<write_type_>::basis_parameters() { return this->basis_; }

template< typename write_type_ >
PetscErrorCode HamiltonianParameters<write_type_>::save_parameters()
{
    common::export_vector_binary(this->prototype_filename(), this->prototype_);
    return Parameters::save_parameters(this->bag_filename);
}
template< typename write_type_ >
PetscErrorCode HamiltonianParameters<write_type_>::init_from_file()
{
    PetscErrorCode e = Parameters::init_from_file(this->bag_filename);
    this->prototype_ = common::import_vector_binary<BasisID>(this->prototype_filename());
    return e;
}
template< typename write_type_ >
PetscErrorCode HamiltonianParameters<write_type_>::register_params()
{
    PetscErrorCode ierr;
    ierr = PetscBagSetName(this->bag,
            "HamiltonianParams",
            "Parameters for finding, and reading the hamiltonian");
    ierr = PetscBagRegisterInt(this->bag, &params->nmax,
            500, "nmax", "Max n value");
    ierr = PetscBagRegisterInt(this->bag, &params->lmax,
            50, "lmax", "Max l value");
    ierr = PetscBagRegisterInt(this->bag, &params->mmax,
            10., "mmax", "Maximum |m| value");
    ierr = PetscBagRegisterString(this->bag, &params->hamiltonian_folder,
            PETSC_MAX_PATH_LEN, "./hamiltonian/",
            "hamiltonian_folder",
            "Where the hamiltonian and the dipole matrix are stored");
    ierr = PetscBagRegisterString(this->bag, &params->basis_bag,
            PETSC_MAX_PATH_LEN, "./basis/BasisParameters.bag",
            "basis_bag_filename",
            "Where the basis parameters bag is");
    return ierr;
}


