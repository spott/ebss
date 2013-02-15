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

class DipoleParameters: public Parameters
{
public:
    DipoleParameters(int argc, const char** argv, MPI_Comm comm): Parameters(comm){
        register_parameters();

        opt.parse(argc, argv);

        get_parameters();

    };

    DipoleParameters(MPI_Comm comm): Parameters(comm){
        //when the parsing will be done elsewhere:
    }

    //call after parsing:
    void get_parameters();

    const std::string dipole_filename() const {return dipole_filename_; };
    const std::string after_filename() const {return after_filename_; };

    PetscReal find_dipole_moment(Mat& dipole, Vec& psi);

    virtual std::string print() const;

    void save_parameters() const;

    PetscReal dt() const { return dt_ / 0.0241888432652; }; //change the fs to A.U.
    PetscReal t_after() const { return t_after_ / 0.0241888432652; }; //ditto

protected:
    ez::ezOptionParser opt;
    void register_parameters();
    bool findDipole;
    std::string dipole_filename_;
    std::string after_filename_;
    PetscReal dt_;
    PetscReal t_after_;
    Vec tmp;

};

PetscReal DipoleParameters::find_dipole_moment(Mat& dipole, Vec& psi)
{
    PetscScalar out;
    if (tmp == PETSC_NULL)
        VecDuplicate(psi, &tmp);
    VecCopy(psi,tmp);
    MatMult(dipole,psi,tmp);
    VecDot(psi,tmp, &out);
    return out.real();
}


void DipoleParameters::get_parameters()
{
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

    if (opt.isSet("-dipole_config"))
    {
        std::string fname;
        opt.get("-dipole_config")->getString(fname);
        if (! opt.importFile(fname.c_str(), '#'))
        {
            std::cout << "file must exist!" << std::endl;
            throw std::exception();
        }
    }
    if (opt.isSet("-dipole_filename"))
    {
        opt.get("-dipole_filename")->getString(dipole_filename_);
        findDipole = true;
    }
    else
        findDipole = false;

    dipole_filename_ = common::absolute_path(dipole_filename_);

    opt.get("-dipole_dt")->getDouble(dt_);
    opt.get("-dipole_t_after")->getDouble(t_after_);
    opt.get("-dipole_after_filename")->getString(after_filename_); 

    tmp = PETSC_NULL;
}

std::string DipoleParameters::print() const
{
    std::ostringstream out;
    out << "dipole_filename: " << dipole_filename_ << std::endl;
    return out.str();
}
void DipoleParameters::save_parameters() const
{
    std::ofstream file;
    file.open(std::string("./Dipole.config"));
    file << "-dipole_filename " << dipole_filename_ << std::endl;
    file.close();
}

void DipoleParameters::register_parameters()
{
    std::string prefix = "-dipole_";
    opt.overview = "Dipole Parameters";
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
            ".1",
            1,
            1,
            0,
            "dt after (in fs)",
            std::string(prefix).append("dt\0").c_str()
           );
    opt.add(
            "10000",
            1,
            1,
            0,
            "t after (in fs)",
            std::string(prefix).append("t_after\0").c_str()
           );
    opt.add(
            "Dipole.dat",
            0,
            1,
            0,
            "dipole filename",
            std::string(prefix).append("filename\0").c_str()
           );
    opt.add(
            "after.dat",
            0,
            1,
            0,
            "after the pulses filename for the dipole",
            std::string(prefix).append("after_filename\0").c_str()
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
