#pragma once

//ebss:
#include<common/common.hpp>
#include<common/math.hpp>
#include<common/parameters/LaserParameters.hpp>

//stl:
#include<sstream>
#include<string>

//petsc:
//#include<petsc.h>

class PulsetrainParameters: public LaserParameters
{
public:
    PulsetrainParameters(int argc, const char** argv, MPI_Comm comm):LaserParameters(argc, argv, comm)
    {
        register_parameters();

        opt.parse(argc, argv);

        get_parameters();
    };

    //call after parsing
    void get_parameters();

    std::vector<double> spacing() const;
    std::vector<double> phase() const;
    virtual PetscScalar efield(PetscReal t);
    PetscReal max_time() const;

    std::string print() const;
protected:
    void register_parameters();
    int num_;
    std::vector<double> spacing_;
    std::vector<double> phase_;
};

std::vector<double> PulsetrainParameters::spacing() const
{
    std::vector<double> s(spacing_.size());

    for(size_t i = 0; i < s.size(); i++)
        s[i] = spacing_[i] / 0.024188843265;

    return s;
}

PetscReal PulsetrainParameters::max_time() const
{
    auto s = spacing();
    double tot;
    for (auto a: s)
        tot += a;

    if (num_ == s.size() + 1)
        return tot + pulse_length();
    else
        return tot > pulse_length() ? tot : pulse_length();
}
std::vector<double> PulsetrainParameters::phase() const { return phase_; }

PetscScalar PulsetrainParameters::efield(PetscReal t)
{
    //we need the AU spacing times:
    auto sp = spacing();
    //so, we use the laser class version, but give it times that correspond to where we are in the pulse.

    PetscScalar ef = 0;
    //First see if the time is below the pulse length:
    double tot_time = 0;
    for (size_t i = 0; i < sp.size(); i++)
    {
        if (t > tot_time)
        {
            cep_ = phase_[i];
            ef += LaserParameters::efield(t - tot_time);
        }
        tot_time += sp[i];
    }
    if (t > tot_time && num_ == sp.size() + 1)
    {
        cep_ = phase_.back();
        ef += LaserParameters::efield(t - tot_time);
    }

    return ef;
}

void PulsetrainParameters::get_parameters()
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

    if (opt.isSet("-pulsetrain_config"))
    {
        std::string fname;
        opt.get("-pulsetrain_config")->getString(fname);
        if (! opt.importFile(fname.c_str(), '#'))
        {
            std::cout << "file must exist!" << std::endl;
            throw std::exception();
        }
    }

    opt.get("-pulsetrain_num")->getInt(num_);
    opt.get("-pulsetrain_spacing")->getDoubles(spacing_);
    opt.get("-pulsetrain_phase")->getDoubles(phase_);
}

std::string PulsetrainParameters::print() const
{

    std::ostringstream out;
    out << LaserParameters::print();
    out << "pulsetrain_num: " << num_ << std::endl;
    out << "pulsetrain_spacing: ";
    for (auto a: spacing_)
        out << a << ",";
    out << std::endl;
    out << "pulsetrain_phase: ";
    for (auto a: phase_)
        out << a << ",";
    out << std::endl;
    return out.str();

}
void PulsetrainParameters::register_parameters()
{
    LaserParameters::register_parameters();

    std::string prefix = "-pulsetrain_";
    opt.overview = "Pulse Train Parameters";
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
            "1",
            0,
            1,
            0,
            "number of pulses",
            std::string(prefix).append("num\0").c_str()
           );
    opt.add(
            "100",
            0,
            1,
            ',',
            "spacing between each pulse (in fs, the last one is repeated)",
            std::string(prefix).append("spacing\0").c_str()
           );
    opt.add(
            "0",
            0,
            1,
            ',',
            "absolute phase of each pulse",
            std::string(prefix).append("phase\0").c_str()
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

