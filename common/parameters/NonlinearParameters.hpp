#pragma once

//ebss:
#include<common/common.hpp>
#include<common/math.hpp>
#include<common/parameters/Parameters.hpp>
#include<common/types.hpp>

//stl:
#include<sstream>
#include<string>

//petsc:
//#include<petsc.h>

class NonlinearParameters: public Parameters
{
public:
    NonlinearParameters(int argc, const char** argv, MPI_Comm comm): Parameters(comm){
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

        if (opt.isSet("-nonlinear_config"))
        {
            std::string fname;
            opt.get("-nonlinear_config")->getString(fname);
            if (! opt.importFile(fname.c_str(), '#'))
            {
                std::cout << "file must exist!" << std::endl;
                throw std::exception();
            }
        }

        if (opt.isSet("-nonlinear_freq_filename"))
        {
            std::string fname;
            opt.get("-nonlinear_freq_filename")->getString(fname);
            if (fname.substr( fname.size() - 3, 3) == "dat")
                freqs_ = common::import_vector_binary<double>(fname);
            else if (fname.substr( fname.size() - 3, 3) == "txt")
                freqs_ = common::import_vector_ascii<double>(fname);
        }
        else
        {
            opt.get("-nonlinear_freq")->getDoubles(freqs_);

            std::vector<double> f;
            opt.get("-nonlinear_energies")->getDoubles(f);
            for( auto a : f)
                freqs_.push_back( units::toEnergy( units::eV( a ) ) );
            f.clear();
            opt.get("-nonlinear_wavelengths")->getDoubles(f);
            for( auto a : f)
                freqs_.push_back( units::toEnergy( units::Meter( a * 1e-9 ) ) );
        }


        opt.get("-nonlinear_img")->getDoubles(imgs_);
        
        std::vector< std::vector<int> > temp;
        opt.get("-nonlinear_chi1")->getMultiInts(temp);
        for( auto a : temp)
        {
            if (a.size() != 1)
                throw("chi1 arguments can only be one element long");
            chi1s_.push_back(a[0]);
        }
        temp.clear();
        opt.get("-nonlinear_chi3")->getMultiInts(temp);
        for( auto a : temp)
        {
            if (a.size() != 3)
                throw("chi3 arguments can only be three element long");
            chi3s_.push_back({{a[0], a[1], a[2]}});
        }
        temp.clear();
        opt.get("-nonlinear_chi5")->getMultiInts(temp);
        for( auto a : temp)
        {
            if (a.size() != 5)
                throw("chi5 arguments can only be three element long");
            chi5s_.push_back({{a[0], a[1], a[2], a[3], a[4]}});
        }

        //sort them:
        std::sort(chi1s_.begin(), chi1s_.end());
        std::sort(chi3s_.begin(), chi3s_.end(), [](std::array<int,3> a,std::array<int,3> b){ 
                if (a[0] < b[0]) return true;
                else if (a[0] >  b[0]) return false;
                else if (a[0] == b[0])
                {
                    if (a[1] < b[1]) return true;
                    else if (a[1] >  b[1]) return false;
                    else if (a[1] == b[1])
                    {
                        if (a[2] < b[2]) return true;
                        else if (a[2] >=  b[2]) return false;
                    }
                }});
        std::sort(chi5s_.begin(), chi5s_.end(), [](std::array<int,5> a,std::array<int,5> b){ 
                if (a[0] < b[0]) return true;
                else if (a[0] >  b[0]) return false;
                else if (a[0] == b[0])
                {
                    if (a[1] < b[1]) return true;
                    else if (a[1] >  b[1]) return false;
                    else if (a[1] == b[1])
                    {
                        if (a[2] < b[2]) return true;
                        else if (a[2] >  b[2]) return false;
                        else if (a[2] == b[2])
                        {
                            if (a[3] < b[3]) return true;
                            else if (a[3] >  b[3]) return false;
                            else if (a[3] == b[3])
                            {
                                if (a[4] < b[4]) return true;
                                else if (a[4] >=  b[4]) return false;
                            }
                        }
                    }
                }});

    };

    std::string print() const;
    void save_parameters() const;

    std::vector<double>& freqs() { return freqs_; };
    std::vector<double>& imgs() {return imgs_; };
    std::vector< int >& chi1s() { return chi1s_; };
    std::vector< std::array<int, 3> >& chi3s() { return chi3s_; };
    std::vector< std::array<int, 5> >& chi5s() { return chi5s_; };

private:
    bool nobound;
    ez::ezOptionParser opt;
    void register_parameters();
    BasisID init_;
    std::string filename_;
    std::vector< int > chi1s_;
    std::vector< std::array<int, 3> > chi3s_;
    std::vector< std::array<int, 5> > chi5s_;
    std::vector< double > freqs_;
    std::vector< double > imgs_;
};

std::string NonlinearParameters::print() const
{
    std::ostringstream out;

    //out << "nonlinear_filename " << filename_ << std::endl;
    for( auto a : chi1s_ )
        out << "nonlinear_chi1 " << a << std::endl;
    for( auto a : chi3s_ )
    {
        out << "nonlinear_chi3 ";
        for (auto b : a)
            out << b << ",";
        out << std::endl;
    }
    for( auto a : chi5s_ )
    {
        out << "nonlinear_chi5 ";
        for (auto b : a)
            out << b << ",";
        out << std::endl;
    }
    out << "nonlinear_img ";
    for (auto a: imgs_)
        out << a << ",";
    out << std::endl;
    out << "nonlinear_freqs ";
    for (auto a : freqs_ )
        out << a << ",";
    out << std::endl;
    return out.str();
}

//TODO fill this out
void NonlinearParameters::save_parameters() const
{
    std::ofstream file;
    file.open(std::string("./Nonlinear.config"));
    file.close();
}

void NonlinearParameters::register_parameters()
{
    std::string prefix = "-nonlinear_";
    opt.overview = "Nonlinear Parameters";
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
            "",
            0,
            1,
            0,
            "the frequencies for a chi1 run (1, 0, -1) are the only allowed values",
            std::string(prefix).append("chi1\0").c_str()
           );
    opt.add(
            "",
            0,
            3,
            ',',
            "the frequencies for a chi3 run (1, 0, -1) are the only allowed values, for example, third harmonic generation is 1,1,1 and n2 will use 1,0,0 or 1,-1,1",
            std::string(prefix).append("chi3\0").c_str()
           );
    opt.add(
            "",
            0,
            5,
            ',',
            "the frequencies for a chi5 run (1, 0, -1) are the only allowed values, for example, third harmonic generation is 1,1,1 and n2 will use 1,0,0 or 1,-1,1",
            std::string(prefix).append("chi5\0").c_str()
           );
    opt.add(
            "",
            0,
            1,
            ',',
            "the frequencies that all \'runs\' will be done at.",
            std::string(prefix).append("freq\0").c_str()
           );
    opt.add(
            ".0001",
            0,
            1,
            ',',
            "the imaginary part of the frequencies",
            std::string(prefix).append("img\0").c_str()
           );
    opt.add(
            "",
            0,
            1,
            ',',
            "the wavelength that all \'runs\' will be done at.",
            std::string(prefix).append("wavelengths\0").c_str()
           );
    opt.add(
            "",
            0,
            1,
            ',',
            "the energies (in au) that all \'runs\' will be done at. This and freq will be added together",
            std::string(prefix).append("energies\0").c_str()
           );
    opt.add(
            "",
            0,
            1,
            0,
            "the filename with the list of frequencies.  If this is used, freq and energies are ignored.",
            std::string(prefix).append("freq_filename\0").c_str()
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
