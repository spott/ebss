#pragma once

#include<vector>
#include<petsc.h>

struct nlm_point {
   int n,l,m;
};

typedef enum { BASIS , DME, PROPAGATE, ETOM } runtype;

typedef struct {
   PetscReal rmax, rmin;
   PetscInt  lmax, nmax, points;
   char basis_folder[PETSC_MAX_PATH_LEN];
} BasisParam;

typedef struct {
   PetscInt nmax, lmax, ammax, basis_size;
   char dme_filename[PETSC_MAX_PATH_LEN];
} DmeParam;

typedef struct {
   PetscReal kmax, kmin, dk;
   char basis_folder[PETSC_MAX_PATH_LEN];
} EtomParam;

typedef struct {
   PetscReal laser_frequency;
   PetscReal laser_efieldstrength;
   PetscReal laser_envelope_size;
   PetscReal laser_cep;
} LaserParam;

typedef struct {
   PetscReal dt;
} PropagateParam;

class parameters{
public:
   //Constructors:
   parameters(MPI_Comm comm, runtype t);
   PetscErrorCode init_from_file(char* filename);

   //save and load the parameters:
   PetscErrorCode save_parameters(const char* filename);
   PetscErrorCode print_parameters();

   //getter methods for calculated values
   const std::vector< nlm_point > vector_prototype();

   //getter methods for stored values:
   const unsigned int nmax() const;
   const unsigned int lmax() const;
   const unsigned int ammax() const;
   const PetscReal rmax() const;
   const PetscReal rmin() const;
   const PetscInt points() const;
   const runtype type() const;
   MPI_Comm comm() const;
   const std::string basis_folder() const;

private:
   PetscErrorCode  register_basis_param(MPI_Comm comm);
   PetscErrorCode  register_dme_param(MPI_Comm comm);
   PetscBag        bag;
   BasisParam      *basis;
   DmeParam        *dme;
   EtomParam       *prop;
   LaserParam      *laser;
   PropagateParam  *etom;
   MPI_Comm        comm_;
   runtype         type_;


};

//// the parameters of the wavefunctions that were generated
//class parameters{

//public:
//parameters(const char *file);
//void out();
//int basis_function_points();
//double rmax();
//double rmin();
//double dr();
//double kmax();
//double kmin();
//double dk();
//double dt();
//double laser_envelope_size();
//double laser_frequency();
//double laser_efieldstrength();
//double laser_cep();
//int nmax();
//int lmax();
//int ammax();
//int nbasis();
//std::string input_folder();
//std::string output_folder();
//std::string filename();
//std::vector< nlm_point >* get_vector_prototype();
//std::vector< nlm_point > make_vector_prototype();
//PetscReal* getBasisFunction(nlm_point state);
//std::string getBasisFunctionFolder();
//std::string getDipoleMatrixFilename();
//std::string getWaveFunctionFilename();

//private:
//double rmax_, rmin_, dr_, kmax_, kmin_, dk_, dt_, laser_intensity_,
//laser_envelope_size_, laser_wavelength_, laser_cep_;
//int nbasis_, lmax_, nmax_, ammax_;  //ammax = max(abs(m))
//std::vector< nlm_point > vector_prototype_;
//std::string output_folder_;
//std::string input_folder_;
//std::string filename_;
//bool folder_exists;
//};

//Constants!
const double PI  = 3.141592653589793238462;
const double HALFPI  = 1.5707963267948966192;
const double SQRT_PI = 1.772453850905516;
const double SQRT_3 = 1.732050807568877;
const double CG10_COEF = 2.046653415892977;
const double SPEED_OF_LIGHT = 137.03599974;
const double NM_TO_AMU = 18.8973;
const double ATTOSECONDS_TO_AMU = 0.0413414;
const double WCM2_TO_AMU = 2.849e-17;
