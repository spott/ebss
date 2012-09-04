#pragma once

//ebss:
#include<common/common.hpp>
#include<common/parameters.hpp>

//Boost:
#include<boost/filesystem.hpp>

//stl
#include<vector>
#include<sstream>
#include<string>

//petsc:
#include<petsc.h>

template<typename T>
class BasisParameters: public Parameters
{
public:
    typedef struct {
        PetscReal rmax, rmin;
        PetscInt  lmax, nmax, points;
        char basis_folder[PETSC_MAX_PATH_LEN];
        char grid_filename[PETSC_MAX_PATH_LEN];
    } basis;

    BasisParameters(MPI_Comm comm): Parameters(comm)
    {
        this->params = new basis;
        p = (void**) params;
        PetscBagCreate(this->comm_, sizeof(basis), &this->bag);
        PetscBagGetData(this->bag, this->p);
        this->register_params();
        PetscBagSetFromOptions(this->bag);

        //see if the folder exists:
        boost::filesystem::create_directories(this->basis_folder());

        this->bag_filename = this->basis_folder().append("/BasisParameters.bag");

        grid_ = new std::vector<T>(this->params->points);
        basis_prototype_ = new std::vector<BasisID>();
    };

    BasisParameters(MPI_Comm comm, std::string filename): Parameters(comm)
    {
        this->bag_filename = filename;
        this->params = new basis;
        p = (void**) params;
        init_from_file();
    };

    //save and init from file.  These might need to be private
    PetscErrorCode save_parameters();
    PetscErrorCode save_parameters(std::string filename);
    PetscErrorCode init_from_file();
    PetscErrorCode init_from_file(std::string filename);

    ~BasisParameters()
    {
        //if it exists or not... write it out.  Either we changed it 
        //(at which point it is important that it gets updated), or we didn't
        //(at which point it doesn't matter).
        this->save_parameters();
        delete grid_;
    };

    //Values!
    PetscReal rmax() const;
    PetscReal rmin() const;
    PetscInt points() const;
    PetscInt nmax() const;
    PetscInt lmax() const;

    //grab the filenames
    std::string basis_folder() const;
    std::string grid_filename() const;
    std::string basis_function_filename(PetscInt n, PetscInt l) const;
    std::string basis_prototype_filename() const;

    //maybe not the best, just handing out references to 
    //the vectors... but, meh, it is just me...
    std::vector<T> * grid();
    std::vector<BasisID> * basis_prototype();
private:
    std::vector<BasisID> *basis_prototype_;
    std::vector<T> *grid_;
    PetscErrorCode register_params();
    basis* params;

};

template <typename T>
PetscErrorCode BasisParameters<T>::init_from_file(std::string filename)
{
    PetscErrorCode e = Parameters::init_from_file(filename);
    grid_ = common::import_vector_binary<T>(this->grid_filename());
    basis_prototype_ = common::import_vector_binary<BasisID>(this->basis_prototype_filename());
    return e;
}

template <typename T>
PetscErrorCode BasisParameters<T>::save_parameters(std::string filename)
{
    std::vector<PetscReal> grid = common::vector_type_change<T, PetscReal>(*this->grid_);
    common::export_vector_binary(this->grid_filename(), &grid);
    common::export_vector_binary(this->basis_prototype_filename(), this->basis_prototype_);
    return Parameters::save_parameters(filename);
}

template <typename T>
PetscErrorCode BasisParameters<T>::init_from_file()
{
    PetscErrorCode e = Parameters::init_from_file(this->bag_filename);
    grid_ = common::import_vector_binary<PetscReal>(this->grid_filename());
    basis_prototype_ = common::import_vector_binary<BasisID>(this->basis_prototype_filename());
    return e;
}

template <typename T>
PetscErrorCode BasisParameters<T>::save_parameters()
{
    std::vector<PetscReal> grid = common::vector_type_change<T, PetscReal>(*this->grid_);
    common::export_vector_binary(this->grid_filename(), &grid);
    common::export_vector_binary(this->basis_prototype_filename(), this->basis_prototype_);
    return Parameters::save_parameters(this->bag_filename);
}

template <typename T>
PetscReal BasisParameters<T>::rmax() const { return params->rmax; }
template <typename T>
PetscReal BasisParameters<T>::rmin() const { return params->rmin; }
template <typename T>
PetscInt BasisParameters<T>::points() const { return params->points; }
template <typename T>
PetscInt BasisParameters<T>::nmax() const { return params->nmax; }
template <typename T>
PetscInt BasisParameters<T>::lmax() const { return params->lmax; }

template <typename T>
std::vector<T> * BasisParameters<T>::grid() { return this->grid_; }
template <typename T>
std::vector<BasisID> * BasisParameters<T>::basis_prototype() { return this->basis_prototype_; }

template <typename T>
std::string BasisParameters<T>::basis_prototype_filename() const
{
    std::stringstream ss;
    std::string s;
    ss << params->basis_folder;
    ss << "prototype.dat";
    ss >> s;
    return s;
}

template <typename T>
std::string BasisParameters<T>::basis_function_filename(PetscInt n, PetscInt l) const
{
    std::stringstream ss;
    std::string s;
    ss << params->basis_folder;
    ss << "/n_" << n << "_l_" << l << ".dat";
    ss >> s;
    return s;
}

template <typename T>
std::string BasisParameters<T>::basis_folder() const 
{  
    std::stringstream ss;
    std::string s;
    ss << params->basis_folder;
    ss >> s;
    return s;
}

template <typename T>
std::string BasisParameters<T>::grid_filename() const
{  
    std::stringstream ss;
    std::string s;
    ss << params->grid_filename;
    ss >> s;
    return s;
}


template <typename T>
PetscErrorCode BasisParameters<T>::register_params()
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
            .000001, "rmin", "Minimum r for grid");
    ierr = PetscBagRegisterInt(this->bag, &params->points,
            10000, "points", "Number of points");
    ierr = PetscBagRegisterString(this->bag, &params->basis_folder,
            PETSC_MAX_PATH_LEN, "./basis/",
            "evectors_folder",
            "Where the evectors and evalues should be stored");
    ierr = PetscBagRegisterString(this->bag, &params->grid_filename,
            PETSC_MAX_PATH_LEN, "./basis/grid.dat",
            "grid array location",
            "Where the grid is");
    return ierr;
}
