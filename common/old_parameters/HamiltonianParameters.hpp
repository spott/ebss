#pragma once

//ebss:
#include<common/common.hpp>
#include<common/parameters/Parameters.hpp>
#include<common/parameters/BasisParameters.hpp>

//stl
#include<vector>
#include<sstream>
#include<string>

//c std lib:
#include<sys/stat.h>

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
        if (flg)
        {
            this->bag_filename = std::string(bagname);
            this->init_from_file();
            this->basis_ = new BasisParameters<write_type_,write_type_>(comm, std::string(this->params->basis_bag));
            //this->print_parameters();
        }
        else
        {
            PetscBagCreate(this->comm_, sizeof(hamiltonian), &this->bag);
            PetscBagGetData(this->bag, (void**)&params);
            this->register_params();
            PetscBagSetFromOptions(this->bag);
            //see if the folder exists:
            mkdir(this->hamiltonian_folder().c_str(), 0775);
            this->bag_filename = this->hamiltonian_folder().append("/Hamiltonian.bag");

            PetscOptionsGetString(PETSC_NULL, "-basis_bag_filename", bagname, PETSC_MAX_PATH_LEN, &flg);
            if (flg)
                this->basis_ = new BasisParameters<write_type_,write_type_>(comm);
            else
            {
                this->basis_ = new BasisParameters<write_type_,write_type_>(comm, std::string(this->params->basis_bag));
            }
        }

        //If we are bigger than the basis we are using, shrink to fit, and print an error message:
        if (this->nmax() > this->basis_->nmax() || this->lmax() > this->basis_->lmax() )
        {
            if (this->rank() == 0) 
            {
                std::cerr << "\e[31m======================ERROR========================" << std::endl;
                std::cerr << "=The basis you are attempting to use is smaller   =" << std::endl;
                std::cerr << "= than you wanted, shrinking to fit.              =" << std::endl;
                std::cerr << "======================ERROR========================\e[0m" << std::endl;
            }
            if (this->nmax() > this->basis_->nmax())
                this->params->nmax = this->basis_->nmax();
            if (this->lmax() > this->basis_->lmax())
                this->params->lmax = this->basis_->lmax();
        }

        //Initialize the prototype file:
        if (this->mmax() == 0 && this->nmax() == this->basis_->nmax() && this->lmax() == this->basis_->lmax())
        {
            this->prototype_.resize(0);
            for(auto i: *(this->basis_->basis_prototype()))
                this->prototype_.push_back(i);
            //std::copy(this->basis_->basis_prototype()->begin(), this->basis_->basis_prototype()->end(), this->prototype_);
            this->prototype_.shrink_to_fit();
        }
        else
        {
            BasisID temp;
            for(int l = 0; l <= this->lmax(); l++)
            {
                for(int m = -( l <= this->mmax() ? l : this->mmax() ); m <= ( l <= this->mmax() ? l : this->mmax() ); m++)
                {
                    for (int n = 1; n <= this->nmax(); n++)
                    {
                        if (n > l)
                        {
                            temp.n = n;
                            temp.l = l;
                            auto el = std::find_if(this->basis_->basis_prototype()->begin(), 
                                                   this->basis_->basis_prototype()->end(),
                                                   [temp](BasisID a)->bool
                                                   { return (a.n == temp.n && a.l == temp.l);});
                            if ( el != this->basis_->basis_prototype()->end())
                                temp.e = el->e;
                            else
                            {
                                std::cout << "The energy wasn't in the basis prototype... I don't know what it is." << std::endl;
                                throw(std::exception());
                            }
                            temp.m = m;
                            this->prototype_.push_back(temp);
                        }
                    }
                }
            }
            this->prototype_.shrink_to_fit();
        }

    }

    //HamiltonianParameters(MPI_Comm comm, BasisParameters<write_type_>* basis): Parameters(comm), basis_(basis)
    //{
        //this->params = new hamiltonian;
        //PetscBagCreate(this->comm_, sizeof(hamiltonian), &this->bag);
        //PetscBagGetData(this->bag, (void**)&params);
        //this->register_params();
        //PetscBagSetFromOptions(this->bag);
    //}

    HamiltonianParameters(MPI_Comm comm, std::string filename): Parameters(comm) 
    {
        this->bag_filename = filename;
        this->params = new hamiltonian;
        this->init_from_file();
        this->basis_ = new BasisParameters<write_type_,write_type_>(comm, std::string(this->params->basis_bag));
    }

    PetscErrorCode save_parameters();
    PetscErrorCode init_from_file();

    ~HamiltonianParameters()
    {
        if (this->rank() == 0) std::cout << "saving parameters" << std::endl;
        this->save_parameters();
        PetscBagDestroy(&this->bag);
    }

    PetscInt nmax() const;
    PetscInt lmax() const;
    PetscInt mmax() const;

    std::string hamiltonian_folder() const;
    std::string dipole_matrix_filename() const;
    void write_dipole_matrix(Mat D);
    Mat read_dipole_matrix();
    std::string energy_eigenvalues_filename() const;
    void write_energy_eigenvalues();
    Vec read_energy_eigenvalues();
    std::string prototype_filename() const;
    BasisParameters<write_type_, write_type_>* basis_parameters();

    std::vector<BasisID> prototype() const;

    PetscErrorCode print_parameters();

private:
    BasisParameters<write_type_, write_type_> *basis_;
    std::vector<BasisID> prototype_;
    hamiltonian* params;
    PetscErrorCode register_params();
};

template<typename write_type_ >
void HamiltonianParameters<write_type_>::write_dipole_matrix(Mat D)
{
    common::petsc_binary_write(this->dipole_matrix_filename(), D, this->comm_);
}

template<typename write_type_ >
Mat HamiltonianParameters<write_type_>::read_dipole_matrix()
{
    return common::petsc_binary_read<Mat>(this->dipole_matrix_filename(), this->comm_);
}

template<typename write_type_ >
void HamiltonianParameters<write_type_>::write_energy_eigenvalues()
{
    Vec ev;
    VecCreate(this->comm_, &ev);
    VecSetType(ev, VECSTANDARD);
    VecSetFromOptions(ev);
    VecSetSizes(ev, PETSC_DECIDE, this->prototype_.size());
    for(size_t i = 0; i < this->prototype_.size(); i++)
    {
        VecSetValue(ev, i, this->prototype_[i].e, INSERT_VALUES);
    }

    VecAssemblyBegin(ev);
    VecAssemblyEnd(ev);

    common::petsc_binary_write(this->energy_eigenvalues_filename(), ev, this->comm_);
}

template<typename write_type_ >
Vec HamiltonianParameters<write_type_>::read_energy_eigenvalues()
{
    return common::petsc_binary_read<Vec>(this->energy_eigenvalues_filename(), this->comm_);
}

template<typename write_type_ >
PetscInt HamiltonianParameters<write_type_>::nmax() const { return this->params->nmax; }

template<typename write_type_ >
PetscInt HamiltonianParameters<write_type_>::lmax() const { return this->params->lmax; }

template<typename write_type_ >
PetscInt HamiltonianParameters<write_type_>::mmax() const { return this->params->mmax; }

template< typename write_type_ >
std::vector<BasisID> HamiltonianParameters<write_type_>::prototype() const { return this->prototype_; }

template< typename write_type_ >
std::string HamiltonianParameters<write_type_>::hamiltonian_folder() const 
{ return std::string(this->params->hamiltonian_folder); } 

template< typename write_type_ >
std::string HamiltonianParameters<write_type_>::dipole_matrix_filename() const
{
    std::stringstream ss;
    std::string s;
    ss << params->hamiltonian_folder;
    ss << "dipole_matrix.dat.gz";
    ss >> s;
    return s;
}


template< typename write_type_ >
std::string HamiltonianParameters<write_type_>::energy_eigenvalues_filename() const
{
    std::stringstream ss;
    std::string s;
    ss << params->hamiltonian_folder;
    ss << "energy_eigenvalues_vector.dat.gz";
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
BasisParameters<write_type_, write_type_>* HamiltonianParameters<write_type_>::basis_parameters() { return this->basis_; }

template< typename write_type_ >
PetscErrorCode HamiltonianParameters<write_type_>::save_parameters()
{
    common::export_vector_binary(this->prototype_filename(), this->prototype_);
    return Parameters::save_parameters(this->bag_filename);
}

template< typename write_type_ >
PetscErrorCode HamiltonianParameters<write_type_>::init_from_file()
{

    PetscErrorCode e = PetscBagCreate(this->comm_, sizeof(hamiltonian), &this->bag);
    PetscBagGetData(this->bag, (void**)&params);
    this->register_params();

    e = Parameters::init_from_file(this->bag_filename);

    this->prototype_ = common::import_vector_binary<BasisID>(this->prototype_filename());
    return e;
}

template< typename write_type_ >
PetscErrorCode HamiltonianParameters<write_type_>::register_params()
{
    PetscErrorCode ierr;
    ierr = PetscBagSetOptionsPrefix(this->bag, "hamiltonian_" );
    ierr = PetscBagSetName(this->bag,
            "HamiltonianParams",
            "Parameters for finding, and reading the hamiltonian");
    ierr = PetscBagRegisterInt(this->bag, &( params->nmax ),
            500, "nmax", "Max n value");
    ierr = PetscBagRegisterInt(this->bag, &( params->lmax ),
            50, "lmax", "Max l value");
    ierr = PetscBagRegisterInt(this->bag, &( params->mmax ),
            0, "mmax", "Maximum |m| value");
    ierr = PetscBagRegisterString(this->bag, &( params->hamiltonian_folder ),
            PETSC_MAX_PATH_LEN, "./hamiltonian/",
            "hamiltonian_folder",
            "Where the hamiltonian and the dipole matrix are stored");
    ierr = PetscBagRegisterString(this->bag, &( params->basis_bag ),
            PETSC_MAX_PATH_LEN, "./basis/BasisParameters.bag",
            "basis_bag_filename",
            "Where the basis parameters bag is");
    return ierr;
}

template< typename write_type_ >
PetscErrorCode HamiltonianParameters<write_type_>::print_parameters()
{
    this->basis_->print_parameters();
    Parameters::print_parameters();
}
