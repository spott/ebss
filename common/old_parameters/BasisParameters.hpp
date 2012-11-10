#pragma once

//ebss:
#include<common/common.hpp>
#include<common/parameters/Parameters.hpp>

//Boost:
//#include<boost/filesystem.hpp>

//stl
#include<vector>
#include<sstream>
#include<string>

//c std lib:
#include<sys/stat.h>

//petsc:
#include<petsc.h>

template<typename compute_type_, typename write_type_ = PetscReal >
class BasisParameters: public Parameters
{
public:

    typedef compute_type_ compute_type;
    typedef write_type_ write_type;
    typedef struct {
        PetscReal rmax, rmin;
        PetscInt  lmax, nmax, points;
        char basis_folder[PETSC_MAX_PATH_LEN];
        char grid_filename[PETSC_MAX_PATH_LEN];
    } basis;

    BasisParameters(MPI_Comm comm): Parameters(comm)
    {
        PetscBool flg = PETSC_FALSE;
        char bagname[PETSC_MAX_PATH_LEN];
        PetscOptionsGetString(PETSC_NULL, "-basis_bag_filename", bagname, PETSC_MAX_PATH_LEN, &flg);

        this->params = new basis;
        if (flg)
        {
            this->bag_filename = std::string(bagname);
            this->init_from_file();
        }
        else
        {
            PetscBagCreate(this->comm_, sizeof(basis), &this->bag);
            PetscBagGetData(this->bag, (void**)&params);

            this->register_params();
            PetscBagSetFromOptions(this->bag);
            //see if the folder exists:
            //std::cout << this->basis_folder() << std::endl;

            mkdir(this->basis_folder().c_str(), 0775);
            //boost::filesystem::create_directories(this->basis_folder());

            this->bag_filename = this->basis_folder().append("/BasisParameters.bag");
            grid_ = std::vector<compute_type>(this->params->points);
            basis_prototype_ = std::vector<BasisID>();
        }
    };

    BasisParameters(MPI_Comm comm, const std::string& filename): Parameters(comm)
    {
        this->bag_filename = filename;
        //this->params = std::shared_ptr<basis>(new basis);
        this->params = new basis;
        this->init_from_file();
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
        std::cerr << "saving parameters" << std::endl;
        this->save_parameters();
        PetscBagDestroy(&this->bag);
        //delete params;
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
    std::string basis_function_filename(BasisID a) const;
    std::string basis_prototype_filename() const;

    //maybe not the best, just handing out references to 
    //the vectors... but, meh, it is just me...
    std::vector<compute_type> * grid();
    std::vector<BasisID> * basis_prototype();
private:
    std::vector<BasisID> basis_prototype_;
    std::vector<compute_type> grid_;
    PetscErrorCode register_params();
    basis* params;

};

template<typename compute_type_, typename write_type_ >
PetscErrorCode BasisParameters<compute_type_, write_type_>::init_from_file(std::string filename)
{
    PetscErrorCode e;
    e = PetscBagCreate(this->comm_, sizeof(basis), &this->bag);
    PetscBagGetData(this->bag, (void**)&params);
    this->register_params();

    Parameters::init_from_file(filename);

    this->grid_ = common::vector_type_change<write_type_, compute_type_>(
            common::import_vector_binary<write_type_>(this->grid_filename())
            );
    this->basis_prototype_ = common::import_vector_binary<BasisID>(this->basis_prototype_filename());
    return e;
}

template<typename compute_type_, typename write_type_ >
PetscErrorCode BasisParameters<compute_type_, write_type_>::init_from_file()
{
    PetscErrorCode e = this->init_from_file(this->bag_filename);
    return e;
}

template<typename compute_type_, typename write_type_ >
PetscErrorCode BasisParameters<compute_type_, write_type_>::save_parameters(std::string filename)
{
    std::vector<write_type_> grid = common::vector_type_change<compute_type_, write_type_>(this->grid_);
    common::export_vector_binary(
            this->grid_filename(), 
            common::vector_type_change<compute_type_, write_type_>(this->grid_)
            );
    common::export_vector_binary(
            this->basis_prototype_filename(),
            this->basis_prototype_);
    return Parameters::save_parameters(filename);
}


template<typename compute_type_, typename write_type_ >
PetscErrorCode BasisParameters<compute_type_, write_type_>::save_parameters()
{
    common::export_vector_binary(
            this->grid_filename(), 
            common::vector_type_change<compute_type_, write_type_>(this->grid_)
            );
    common::export_vector_binary(
            this->basis_prototype_filename(), 
            this->basis_prototype_);
    return Parameters::save_parameters(this->bag_filename);
}

template<typename compute_type_, typename write_type_ >
PetscReal BasisParameters<compute_type_, write_type_>::rmax() const { return params->rmax; }
template<typename compute_type_, typename write_type_ >
PetscReal BasisParameters<compute_type_, write_type_>::rmin() const { return params->rmin; }
template<typename compute_type_, typename write_type_ >
PetscInt BasisParameters<compute_type_, write_type_>::points() const { return params->points; }
template<typename compute_type_, typename write_type_ >
PetscInt BasisParameters<compute_type_, write_type_>::nmax() const { return params->nmax; }
template<typename compute_type_, typename write_type_ >
PetscInt BasisParameters<compute_type_, write_type_>::lmax() const { return params->lmax; }

template<typename compute_type_, typename write_type_ >
std::vector<compute_type_> * BasisParameters<compute_type_, write_type_>::grid() { return &this->grid_; }
template<typename compute_type_, typename write_type_ >
std::vector<BasisID> * BasisParameters<compute_type_, write_type_>::basis_prototype() { return &this->basis_prototype_; }

template<typename compute_type_, typename write_type_ >
std::string BasisParameters<compute_type_, write_type_>::basis_prototype_filename() const
{
    std::stringstream ss;
    std::string s;
    ss << params->basis_folder;
    ss << "prototype.dat";
    ss >> s;
    return s;
}

template<typename compute_type_, typename write_type_ >
std::string BasisParameters<compute_type_, write_type_>::basis_function_filename(BasisID a) const
{
    std::stringstream ss;
    std::string s;
    ss << params->basis_folder;
    ss << "/n_" << a.n << "_l_" << a.l << ".dat";
    ss >> s;
    return s;
}
template<typename compute_type_, typename write_type_ >
std::string BasisParameters<compute_type_, write_type_>::basis_function_filename(PetscInt n, PetscInt l) const
{
    std::stringstream ss;
    std::string s;
    ss << params->basis_folder;
    ss << "/n_" << n << "_l_" << l << ".dat";
    ss >> s;
    return s;
}

template<typename compute_type_, typename write_type_ >
std::string BasisParameters<compute_type_, write_type_>::basis_folder() const 
{  
    std::stringstream ss;
    std::string s;
    ss << params->basis_folder;
    ss >> s;
    return s;
}

template<typename compute_type_, typename write_type_ >
std::string BasisParameters<compute_type_, write_type_>::grid_filename() const
{
    std::stringstream ss;
    std::string s;
    ss << params->grid_filename;
    ss >> s;
    return s;
}


template<typename compute_type_, typename write_type_ >
PetscErrorCode BasisParameters<compute_type_, write_type_>::register_params()
{
    PetscErrorCode ierr;
    ierr = PetscBagSetOptionsPrefix(this->bag, "basis_" );
    ierr = PetscBagSetName(bag,
            "BasisParams",
            "Parameters for finding, and reading the basis state");
    ierr = PetscBagRegisterInt(bag, &(params->nmax),
            500, "nmax", "Max n value");
    ierr = PetscBagRegisterInt(this->bag, &(params->lmax),
            50, "lmax", "Max l value");
    ierr = PetscBagRegisterReal(this->bag, &(params->rmax),
            1000., "rmax", "Maximum r for grid");
    ierr = PetscBagRegisterReal(this->bag, &(params->rmin),
            .000001, "rmin", "Minimum r for grid");
    ierr = PetscBagRegisterInt(this->bag, &(params->points),
            10000, "points", "Number of points");
    ierr = PetscBagRegisterString(this->bag, &(params->basis_folder),
            PETSC_MAX_PATH_LEN, "./basis/",
            "evectors_folder",
            "Where the evectors and evalues should be stored");
    ierr = PetscBagRegisterString(this->bag, &(params->grid_filename),
            PETSC_MAX_PATH_LEN, "./basis/grid.dat",
            "grid_location",
            "Where the grid is");
    return ierr;
}
