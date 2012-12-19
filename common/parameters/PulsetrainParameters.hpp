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
        spacing();
    };


    //call after parsing
    void get_parameters();

    bool in_pulse(PetscReal t);
    std::vector<double> spacing();
    std::vector<double> phase() const;
    virtual PetscScalar efield(PetscReal t);
    PetscReal max_time();
    PetscReal next_pulse_start(PetscReal t);

    void save_parameters() const;
    std::string print() const;
protected:
    void register_parameters();
    std::vector<double> spacing_;
    std::vector<double> phase_;
    std::vector<double> s_au;  //spacing in au
};

std::vector<double> PulsetrainParameters::spacing()
{
    if (s_au.size() == spacing_.size())
        return s_au;


    s_au.resize(spacing_.size());
    for(size_t i = 0; i < s_au.size(); i++)
        s_au[i] = spacing_[i] / 0.024188843265; //translate to AU

    return s_au;
}

bool PulsetrainParameters::in_pulse(PetscReal t)
{
    double tot_time = 0;
    for (size_t i = 0; i < s_au.size(); i++)
    {
        if ( t < tot_time + s_au[i] && t > pulse_length() + tot_time)
            return false;
        tot_time += s_au[i];
    }
    return true;
}

PetscReal PulsetrainParameters::next_pulse_start(PetscReal t)
{
    PetscReal start = 0;
    PetscReal total = 0;
    for (size_t i = 0; i < s_au.size(); i++)
    {
        total += s_au[i];
        if (t < total)
            return total;
    }
    return 0;
}


PetscReal PulsetrainParameters::max_time()
{
    double tot;
    for (auto a: s_au)
        tot += a; //one pulse per, the "spacing" is peak to peak.

    return tot + pulse_length();
}

std::vector<double> PulsetrainParameters::phase() const { return phase_; }

PetscScalar PulsetrainParameters::efield(PetscReal t)
{
    //find the start of the pulse:
    PetscReal start = 0;
    PetscReal total = 0;
    for (size_t i = 0; i < s_au.size(); i++)
    {
        if (t > total)
            start = total;
        total += s_au[i];
    }
    if (t > total)
        start = total;

    if (in_pulse(t))
    {
        PetscScalar efield = std::sqrt(this->intensity());
        efield *= envelope(t, start) * std::sin( this->frequency() * t);
        return efield;
    }
    else
        return 0;
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

    opt.get("-pulsetrain_spacing")->getDoubles(spacing_);
    opt.get("-pulsetrain_phase")->getDoubles(phase_);
}

void PulsetrainParameters::save_parameters() const
{
    std::ofstream file;
    file.open(std::string("./Pulsetrain.config"));
    file << "-pulsetrain_spacing ";
    for (size_t i = 0; i < spacing_.size(); i++)
    {
        if (i == 0)
            file << spacing_[i];
        else
            file << "," << spacing_[i];
    }
    file << std::endl;
    file << "-pulsetrain_phase ";
    for (size_t i = 0; i < phase_.size(); i++)
    {
        if (i == 0)
            file << phase_[i];
        else
            file << "," << phase_[i];
    }
    file << std::endl;
    file.close();
}

std::string PulsetrainParameters::print() const
{
    std::ostringstream out;
    out << LaserParameters::print();
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

