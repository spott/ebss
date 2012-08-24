#include<iomanip>
#include<petsc.h>
#include<vector>
#include<common/parameters.hpp>


parameters::parameters(MPI_Comm comm, runtype t)
{
   this->type_ = t;
   this->comm_ = comm;
   if ( t == BASIS )
   {
      PetscBagCreate(comm, sizeof(BasisParam), &this->bag);
      PetscBagGetData(this->bag, (void **) &basis);
      this->register_basis_param(comm);
      dme = NULL;
      prop = NULL;
      laser = NULL;
      etom = NULL;
   }
   else if ( t == DME )
   {
      PetscBagCreate(comm, sizeof(DmeParam), &this->bag);
      PetscBagGetData(this->bag, (void **) &dme);
      this->register_dme_param(comm);
      basis = NULL;
      prop = NULL;
      laser = NULL;
      etom = NULL;
   }
   PetscBagSetFromOptions(this->bag);
}

PetscErrorCode parameters::init_from_file(char* filename)
{
   PetscViewer viewer;
   PetscViewerBinaryOpen(this->comm_,filename,FILE_MODE_READ,&viewer);
   PetscBagLoad(viewer, &this->bag);
   if ( this->type_ == BASIS )
   {
      PetscBagGetData(this->bag, (void **) &basis);
   }
   else if ( this->type_ == DME )
   {
      PetscBagGetData(this->bag, (void **) &dme);
   }
   PetscBagSetFromOptions(this->bag);
   PetscViewerDestroy(&viewer);
   return 0;
}

PetscErrorCode parameters::register_dme_param(MPI_Comm comm)
{
   PetscErrorCode ierr;
   ierr = PetscBagSetName(this->bag,
                          "DmeParams",
                          "Parameters for finding the dme state");
   ierr = PetscBagRegisterInt(this->bag, &dme->nmax,
                              500, "nmax", "Max n value");
   ierr = PetscBagRegisterInt(this->bag, &dme->lmax,
                              50, "lmax", "Max l value");
   ierr = PetscBagRegisterInt(this->bag, &dme->ammax,
                              0, "ammax", "Max m value");
   ierr = PetscBagRegisterString(this->bag, &dme->dme_filename,
                                 PETSC_MAX_PATH_LEN, "./dme.dat","dme_filename",
                                 "the base filename for the dme");
   return ierr;
}

PetscErrorCode parameters::register_basis_param(MPI_Comm comm)
{
   PetscErrorCode ierr;
   ierr = PetscBagSetName(this->bag,
                          "BasisParams",
                          "Parameters for finding the basis state");
   ierr = PetscBagRegisterInt(this->bag, &basis->nmax,
                              500, "nmax", "Max n value");
   ierr = PetscBagRegisterInt(this->bag, &basis->lmax,
                              50, "lmax", "Max l value");
   ierr = PetscBagRegisterReal(this->bag, &basis->rmax,
                               1000., "rmax", "Maximum r for grid");
   ierr = PetscBagRegisterReal(this->bag, &basis->rmin,
                               .1, "rmin", "Minimum r for grid");
   ierr = PetscBagRegisterInt(this->bag, &basis->points,
                               10000, "points", "Number of points");
   ierr = PetscBagRegisterString(this->bag, &basis->basis_folder,
                                 PETSC_MAX_PATH_LEN, "./basis/",
                                 "evectors_folder",
                          "Where the evectors and evalues should be stored");

   return ierr;
}

PetscErrorCode parameters::save_parameters(const char* filename)
{
   PetscErrorCode ierr;
   PetscViewer viewer;
   ierr = PetscViewerBinaryOpen(this->comm_,filename,FILE_MODE_WRITE,&viewer);
   ierr = PetscBagView(this->bag, viewer);
   ierr = PetscViewerDestroy(&viewer);

   return ierr;
}

PetscErrorCode parameters::print_parameters()
{
   PetscErrorCode ierr;
   ierr = PetscBagView(this->bag,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
   return (0);
}

const unsigned int parameters::nmax() const
{
   if (this->type_ == DME)
      return this->dme->nmax;
   else if (this->type_ == BASIS)
      return this->basis->lmax;
   else
      throw(std::exception());
}

const unsigned int parameters::lmax() const
{
   if (this->type_ == DME)
      return this->dme->lmax;
   else if (this->type_ == BASIS)
      return this->basis->lmax;
   else
      throw(std::exception());
}

const unsigned int parameters::ammax() const
{
   if (this->type_ == DME)
      return this->dme->ammax;
   else
      throw(std::exception());
}

const PetscReal parameters::rmax() const
{
   if (this->type_ == BASIS)
      return this->basis->rmax;
   else
      throw(std::exception());
}
const PetscReal parameters::rmin() const
{
   if (this->type_ == BASIS)
      return this->basis->rmin;
   else
      throw(std::exception());
}
const PetscInt parameters::points() const
{
   if (this->type_ == BASIS)
      return this->basis->points;
   else
      throw(std::exception());
}

const runtype parameters::type() const
{
   return this->type_;
}

MPI_Comm parameters::comm() const
{
   return this->comm_;
}

const std::string parameters::basis_folder() const
{
   if (this->type_ == BASIS)
      return std::string(this->basis->basis_folder);
   else
      throw(std::exception());
}

const std::vector< nlm_point > parameters::vector_prototype()
{
   if (this->type_ != DME || this->type_ != BASIS)
      throw(std::exception());

   int lmax = (this->type_ == BASIS) ? this->basis->lmax : this->dme->lmax;
   int nmax = (this->type_ == BASIS) ? this->basis->nmax : this->dme->nmax;
   std::vector< nlm_point > vec_prototype(0);
   if (this->type_ == DME)
      this->dme->basis_size = 0;
   for (int l = 0; l <= lmax; l++)
   {
      for (int n = l + 1; n <= nmax; n++)
      {
         nlm_point tmp;
         tmp.n = n;
         tmp.l = l;
         tmp.m = 0;
         vec_prototype.push_back( tmp );
         if (this->type_ == DME)
            this->dme->basis_size++;
      }
   }
   return vec_prototype;
}
//const std::vector< nlm_point > parametersvector_prototype(int nmax, int lmax)
//{
//std::vector< nlm_point > vec_prototype(0);
//for (int l = 0; l <= lmax; l++)
//{
//for (int n = l + 1; n <= nmax; n++)
//{
//nlm_point tmp;
//tmp.n = n;
//tmp.l = l;
//tmp.m = m;
//vec_prototype.push_back( tmp );
//}
//}
//return vec_prototype;
//}
//parameters::parameters(const char *file)
//{
//this->rmax_ = 0;     this->rmin_ = 0;   this->dr_ = 0;     this->kmax_ = 0;
//this->kmin_ = 0;     this->dk_ = 0;     this->dt_ = 0;     this->laser_intensity_ = 0;
//this->laser_envelope_size_ = 0; this->laser_wavelength_ = 0; this->laser_cep_ = 0;
//this->nbasis_ = 0;   this->lmax_ = 0;   this->nmax_ = 0;   this->ammax_ = 0;
//this->folder_exists = false;
//this->filename_ = std::string(file);
////std::cout << this->folder <<std::endl;
//std::ifstream ifile(this->filename_.c_str(),std::ios::in);
//if(ifile.is_open())
//{
//std::string tmp;
//while(!ifile.eof())
//{
//ifile >> tmp;
//if (tmp == "rmax")
//ifile >> this->rmax_;
//else if (tmp == "rmin")
//ifile >> this->rmin_;
//else if (tmp == "dr")
//ifile >> this->dr_;
//else if (tmp == "kmax")
//ifile >> this->kmax_;
//else if (tmp == "kmin")
//ifile >> this->kmin_;
//else if (tmp == "dk")
//ifile >> this->dk_;
//else if (tmp == "nmax")
//ifile >> this->nmax_;
//else if (tmp == "lmax")
//ifile >> this->lmax_;
//else if (tmp == "ammax")
//ifile >> this->ammax_;
//else if (tmp == "dt")
//ifile >> this->dt_;
//else if (tmp == "envelope_size")
//ifile >> this->laser_envelope_size_;
//else if (tmp == "laser_wavelength")
//ifile >> this->laser_wavelength_;
//else if (tmp == "laser_intensity")
//ifile >> this->laser_intensity_;
//else if (tmp == "laser_cep")
//ifile >> this->laser_cep_;
//else if (tmp == "output_folder")
//ifile >> this->output_folder_;
//else if (tmp == "input_folder")
//ifile >> this->input_folder_;
//else
//std::cerr << "parameter value not supported: |" << tmp << "|" << std::endl;
//};

//ifile.close();
//}
//else
//{
//std::cout << "parameters file isn't what we expect:" << std::endl;
//std::cout << this->filename_ << std::endl;
//throw(std::exception());
//}
//}

//PetscReal* parameters::getBasisFunction(nlm_point state)
//{
//if (state.n == 0)
//std::cerr << "This function is asking for n = 0, why?" << std::endl;
//int npoints = basis_function_points();
//PetscReal* psi = new PetscReal[npoints];

//std::stringstream file(std::ios::out);
//file << getBasisFunctionFolder() << "/evectors_n" << std::setw(4) << std::setfill('0') << state.n << "_l" << std::setw(3) << std::setfill('0') << state.l;
//std::ifstream ifile(file.str().c_str(), std::ios::binary | std::ios::in);
//if(ifile.is_open())
//{
//ifile.read((char *)(psi), sizeof(PetscReal)*npoints);
//ifile.close();
//}
//else
//{
//std::cerr << "file didn't open !" << std::endl;
//std::cerr << file.str() << std::endl;
//}

////check info:
//int errors=0;
//for (int i = 0; i < npoints; i++)
//{
//if ( !(psi[i] <= DBL_MAX && psi[i] >= -DBL_MAX) )
//{
//std::cerr << "the wavefunction is wrong" << std::endl;
//errors++;
//}
//}
//if (errors != 0 )
//std::cerr << "file: " << file.str() << " has " << errors << " errors in it" << std::endl;

//return psi;
//}

//std::string parameters::getDipoleMatrixFilename()
//{
//std::stringstream folder(this->input_folder_);
//folder << this->input_folder_;
//folder << "dme_n" << this->nmax_ << "_l" << this->lmax_ << "_m" << this->ammax_ << ".dat";
//return folder.str();
//}
//std::string parameters::getWaveFunctionFilename()
//{
//std::stringstream folder(this->output_folder_);
//folder << this->output_folder_;
//folder << "wf_final_n" << this->nmax_ << "_l" << this->lmax_ << "_m" << this->ammax_ << "_dt" << this->dt_ << ".dat";
//return folder.str();
//}

//std::vector< nlm_point >* parameters::get_vector_prototype()
//{
//if (this->vector_prototype_.empty())
//return &this->vector_prototype_;
//this->make_vector_prototype();
//return &this->vector_prototype_;
//}

//std::vector<nlm_point> parameters::make_vector_prototype()
//{
//std::vector< nlm_point > vec_prototype(0);
//this->nbasis_ = 0;
//for (int l = 0; l <= this->lmax_; l++)
//{
//for (int n = l + 1; n <= this->nmax_; n++)
//{
//for (int m = (this->ammax_ < l ? -this->ammax_ : -l ) ; m <= (this->ammax_ < l ? this->ammax_ : l) ; m++)
//{
//nlm_point tmp;
//tmp.n = n;
//tmp.l = l;
//tmp.m = m;
//vec_prototype.push_back( tmp );
//this->nbasis_++;
//}
//}
//}
//return vec_prototype;
//}

//std::string parameters::getBasisFunctionFolder()
//{
//std::stringstream folder(this->input_folder_);
//folder << this->input_folder_;
//folder << "/basis/";

//return folder.str();
//}

//double parameters::rmax() {return this->rmax_;}
//double parameters::rmin() {return this->rmin_;}
//double parameters::dr() {return this->dr_;}
//double parameters::kmax() {return this->kmax_;}
//double parameters::kmin() {return this->kmin_;}
//double parameters::dk() {return this->dk_;}
//int parameters::nmax() {return this->nmax_;}
//int parameters::nbasis() {this->make_vector_prototype(); return this->nbasis_;}
//int parameters::lmax() {return this->lmax_;}
//int parameters::ammax() {return this->ammax_;}
//std::string parameters::input_folder() {return this->input_folder_;}
//std::string parameters::output_folder() {return this->output_folder_;}
//std::string parameters::filename() {return this->filename_;}

////Convert from real values (SI units, kind of) to Atomic Units:

////time in atoseconds to time in amu:
//double parameters::dt()
//{
////return this->dt_ * ATTOSECONDS_TO_AMU;
//return this->dt_;
//}

////number of cycles to amu
//double parameters::laser_envelope_size()
//{
//return this->laser_frequency()/ this->laser_envelope_size_ ;
//}

//// laser wavelength in nm to angular frequency in amu
//double parameters::laser_frequency()
//{
////double lambda = this->laser_wavelength_ * NM_TO_AMU;
////return SPEED_OF_LIGHT/ (lambda * 2 * PI);
//return 45.5896 / this->laser_wavelength_;
//}
//double parameters::laser_efieldstrength()
//{
////return std::sqrt((this->laser_intensity_ * WCM2_TO_AMU * 8. * PI)/SPEED_OF_LIGHT) ;
//return std::sqrt(this->laser_intensity_ / 3.5101e+16 ) ;
//}

//double parameters::laser_cep() {return this->laser_cep_;}

//int parameters::basis_function_points()
//{
//return int( ( this->rmax_ - this->rmin_)/ this->dr_ ) + 1;
//}

//void parameters::out()
//{
//std::cout << "rmax:   " << this->rmax_ << std::endl;
//std::cout << "rmin:   " << this->rmin_ << std::endl;
//std::cout << "dr:     " << this->dr_ << std::endl;
//std::cout << "kmax:   " << this->kmax_ << std::endl;
//std::cout << "kmin:   " << this->kmin_ << std::endl;
//std::cout << "dk:     " << this->dk_ << std::endl;
//std::cout << "lmax:   " << this->lmax_ << std::endl;
//std::cout << "nmax:   " << this->nmax_ << std::endl;
//std::cout << "|mmax|: " << this->ammax_ << std::endl;
//std::cout << "dt:     " << this->dt_ << std::endl;
//std::cout << "output_folder: " << this->output_folder_ << std::endl;
//std::cout << "input_folder: " << this->input_folder_ << std::endl;
//std::cout << "filename: " << this->filename_ << std::endl;
//}
