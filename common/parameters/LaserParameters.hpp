#pragma once

//ebss:
#include<common/common.hpp>
#include<common/math.hpp>
#include<common/parameters/Parameters.hpp>

//stl:
#include<sstream>
#include<string>

//petsc:
//#include<petsc.h>

class LaserParameters: public Parameters
{
public:
    LaserParameters(int argc, const char** argv, MPI_Comm comm): Parameters(comm){
        register_parameters();

        opt.parse(argc, argv);

        if (opt.isSet("+d")) {
            std::string pretty;
            opt.prettyPrint(pretty);
            std::cout << pretty;
        }
        if (opt.isSet("-h")) {
            std::string usage;
            opt.getUsage(usage,80,ez::ezOptionParser::ALIGN);
            std::cout << usage;
        }

        if (opt.isSet("-laser_config"))
        {
            std::string fname;
            opt.get("-laser_config")->getString(fname);
            if (! opt.importFile(fname.c_str(), '#'))
            {
                std::cout << "file must exist!" << std::endl;
                throw std::exception();
            }
        }

        opt.get("-laser_lambda")->getDouble(lambda_);
        opt.get("-laser_intensity")->getDouble(intensity_);
        opt.get("-laser_cep")->getDouble(cep_);
        opt.get("-laser_cycles")->getDouble(cycles_);
        opt.get("-laser_dt")->getDouble(dt_);
        opt.get("-laser_dt_after")->getDouble(dt_after_);
        opt.get("-laser_t_after")->getDouble(t_after_);
        opt.get("-laser_laser_filename")->getString(laser_filename_);
    };

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

    std::string print() const;

    void save_parameters() const;

private:
    ez::ezOptionParser opt;
    void register_parameters();
    double lambda_;
    double intensity_;
    double cep_;
    double cycles_;
    double dt_;
    double dt_after_;
    double t_after_;
    std::string laser_filename_;

};

std::string LaserParameters::print() const
{
    std::ostringstream out;
    out << "laser_lambda: " << lambda_ << std::endl;
    out << "laser_intensity: " << intensity_ << std::endl;
    out << "laser_cep: " << cep_ << std::endl;
    out << "laser_cycles" << cycles_ << std::endl;
    out << "laser_dt" << dt_ << std::endl;
    out << "laser_dt_after" << dt_after_ << std::endl;
    out << "laser_t_after" << t_after_ << std::endl;
    out << "laser_filename" << laser_filename_ << std::endl;
    return out.str();
}
void LaserParameters::save_parameters() const
{
    std::ofstream file;
    file.open(std::string("./Laser.config"));
    file << "-laser_lambda " << lambda_ << std::endl;
    file << "-laser_intensity " << intensity_ << std::endl;
    file << "-laser_cep " << cep_ << std::endl;
    file << "-laser_cycles " << cycles_ << std::endl;
    file << "-laser_dt " << dt_ << std::endl;
    file << "-laser_dt_after " << dt_after_ << std::endl;
    file << "-laser_t_after " << t_after_ << std::endl;
    file << "-laser_laser_filename " << laser_filename_ << std::endl;
    file.close();
}

PetscReal LaserParameters::lambda() const { return (lambda_ / 5.29177206e-2); }
PetscReal LaserParameters::frequency() const { return (math::C * 2 * math::PI) / (this->lambda() ); }
PetscReal LaserParameters::intensity() const
{ return intensity_/3.5094452e16; }
PetscReal LaserParameters::cep() const { return cep_; }
PetscReal LaserParameters::cycles() const { return cycles_; }
PetscReal LaserParameters::dt() const { return dt_; }
PetscReal LaserParameters::dt_after() const { return dt_after_; }
PetscReal LaserParameters::t_after() const { return t_after_; }
std::string LaserParameters::laser_filename() const { return std::string(laser_filename_); }

PetscScalar 
LaserParameters::efield(PetscReal t)
{
    if (t * this->frequency()/this->cycles() > math::PI)
        return 0.0;
    PetscReal efield = std::sqrt(this->intensity());
    return efield 
        * std::pow( std::sin( this->frequency() * t / (this->cycles() * 2) ) ,2) 
        * std::sin( this->frequency() * t + this->cep() );
}


void LaserParameters::register_parameters()
{
    std::string prefix = "-laser_";
    opt.overview = "Laser Parameters";
    opt.add(
            "", // Default.
            0, // Required?
            0, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "Display usage instructions.", // Help description.
            "-h",     // Flag token. 
            "-help",  // Flag token.
            "--help", // Flag token.
            "--usage" // Flag token.
           );
    opt.add(
            "800",
            0,
            1,
            0,
            "wavelength of the laser in nm",
            std::string(prefix).append("lambda\0").c_str()
           );
    opt.add(
            "10e12",
            0,
            1,
            0,
            "intensity of laser in W/cm^2",
            std::string(prefix).append("intensity\0").c_str()
           );
    opt.add(
            "0",
            0,
            1,
            0,
            "Carrier Envelope Phase",
            std::string(prefix).append("cep\0").c_str()
           );
    opt.add(
            "10",
            0,
            1,
            0,
            "number of cycles of the carrier frequency",
            std::string(prefix).append("cycles\0").c_str()
           );
    opt.add(
            "0.01",
            0,
            1,
            0,
            "timestep durring the pulse (in A.U.)",
            std::string(prefix).append("dt\0").c_str()
           );
    opt.add(
            "0.01",
            0,
            1,
            0,
            "timestep after the pulse (in A.U.)",
            std::string(prefix).append("dt_after\0").c_str()
           );
    opt.add(
            "0",
            0,
            1,
            0,
            "time after the pulse (in A.U.)",
            std::string(prefix).append("t_after\0").c_str()
           );
    opt.add(
            "./efield.dat",
            0,
            1,
            0,
            "filename to put the laser",
            std::string(prefix).append("laser_filename\0").c_str()
           );
    opt.add(
            "",
            0,
            1,
            0,
            "Config file to import",
            std::string(prefix).append("config\0").c_str()
           );
    opt.add(
            "", // Default.
            0, // Required?
            0, // Number of args expected.
            0, // Delimiter if expecting multiple args.
            "Print all inputs and categories for debugging.", // Help description.
            "+d",
            "--debug"     // Flag token. 
           );
}
