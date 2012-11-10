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

class LaserParameters: public Parameters
{
public:
    typedef struct {
        PetscReal lambda,
                  intensity,
                  cep,
                  cycles;
        PetscReal dt,
                  dt_after,
                  t_after;
        char laser_filename[PETSC_MAX_PATH_LEN];
    } laser;

    LaserParameters(MPI_Comm comm): Parameters(comm)
    {
        PetscBool flg = PETSC_FALSE;
        char bagname[PETSC_MAX_PATH_LEN];
        PetscOptionsGetString(PETSC_NULL, "-laser_bag_filename", bagname, PETSC_MAX_PATH_LEN, &flg);

        this->params = new laser;
        //p = (void**) params;
        if (flg)
        {
            this->bag_filename = std::string(bagname);
            this->init_from_file();
        }
        else
        {
            PetscBagCreate(this->comm_, sizeof(laser), &this->bag);
            PetscBagGetData(this->bag, (void**)&params);
            this->register_params();
            PetscBagSetFromOptions(this->bag);
            //see if the folder exists:
            this->bag_filename = std::string("./Laser.bag");
        }
    }

    LaserParameters(MPI_Comm comm, std::string filename): Parameters(comm)
    {
    }

    PetscReal lambda() const;
    PetscReal frequency() const;
    PetscReal intensity() const;
    PetscReal cep() const;
    PetscReal cycles() const;
    PetscReal dt() const;
    PetscReal dt_after() const;
    PetscReal t_after() const;
    std::string laser_filename() const;
    PetscScalar efield(PetscReal t);

    PetscErrorCode init_from_file();
    PetscErrorCode save_parameters();


private:
    laser* params;
    PetscErrorCode register_params();

};
PetscErrorCode LaserParameters::init_from_file()
{
    PetscErrorCode e = PetscBagCreate(this->comm_, sizeof(laser), &this->bag);
    e = PetscBagGetData(this->bag, (void**)&params);
    this->register_params();

    e = Parameters::init_from_file(this->bag_filename);
    return e;
}

PetscErrorCode LaserParameters::save_parameters()
{
    return Parameters::save_parameters(this->bag_filename);
}

PetscReal LaserParameters::lambda() const { return (this->params->lambda / 5.29177206e-2); }
PetscReal LaserParameters::frequency() const { return (math::C * 2 * math::PI) / (this->lambda() ); }
PetscReal LaserParameters::intensity() const
{ return this->params->intensity/3.5094452e16; }
PetscReal LaserParameters::cep() const { return this->params->cep; }
PetscReal LaserParameters::cycles() const { return this->params->cycles; }
PetscReal LaserParameters::dt() const { return this->params->dt; }
PetscReal LaserParameters::dt_after() const { return this->params->dt_after; }
PetscReal LaserParameters::t_after() const { return this->params->t_after; }
std::string LaserParameters::laser_filename() const { return std::string(this->params->laser_filename); }

PetscScalar 
LaserParameters::efield(PetscReal t)
{
    if (t * this->frequency()/this->cycles() > math::PI)
        return 0.0;
    PetscReal efield = std::sqrt(this->intensity());
    return efield 
        * std::pow( std::sin( this->frequency() * t / this->cycles() ) ,2) 
        * std::sin( this->frequency() * t + this->cep() );
}

PetscErrorCode 
LaserParameters::register_params()
{
    PetscErrorCode ierr;
    ierr = PetscBagSetOptionsPrefix(this->bag, "laser_" );
    ierr = PetscBagSetName(this->bag,
            "LaserParams",
            "Parameters for the laser");
    ierr = PetscBagRegisterReal(this->bag, &params->lambda,
            800, "lambda", "wavelength of the laser in nm");
    ierr = PetscBagRegisterReal(this->bag, &params->intensity,
            10e12, "intensity", "intensity of the laser in W/cm^2");
    ierr = PetscBagRegisterReal(this->bag, &params->cep,
            0, "cep", "carrier envelope phase");
    ierr = PetscBagRegisterReal(this->bag, &params->cycles,
            10, "cycles", "number of cycles in the pulse");
    ierr = PetscBagRegisterReal(this->bag, &params->dt,
            0.01, "dt", "timestep durring the pulse");
    ierr = PetscBagRegisterReal(this->bag, &params->dt_after,
            0.01, "dt_after", "the timestep after the pulse");
    ierr = PetscBagRegisterReal(this->bag, &params->t_after,
            0, "t_after", "the time after the pulse");
    ierr = PetscBagRegisterString(this->bag, &params->laser_filename,
            PETSC_MAX_PATH_LEN, "./laser.dat",
            "laser_filename",
            "The file to put the laser field in.  If it already exists, the file the laser field is read from (the rest of the parameters are ignored)");
    return ierr;
}
