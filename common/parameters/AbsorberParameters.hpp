#pragma once

#include<common/common.hpp>
#include<common/parameters/Parameters.hpp>

class AbsorberParameters: public Parameters
{
public:
    AbsorberParameters(int argc, const char** argv, MPI_Comm comm): Parameters(comm){
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

        if (opt.isSet("-absorber_config"))
        {
            std::string fname;
            opt.get("-absorber_config")->getString(fname);
            if (! opt.importFile(fname.c_str(), '#'))
            {
                std::cout << "file must exist!" << std::endl;
                throw std::exception();
            }
        }

        opt.get("-absorber_n_size")->getInt(n_size_);
        opt.get("-absorber_l_size")->getInt(l_size_);
        opt.get("-absorber_m_size")->getInt(m_size_);
        opt.get("-absorber_cos_factor")->getDouble(cos_factor_);
    };

    //AbsorberParameters(std::string filename, MPI_Comm comm): Parameters(comm) {
        //init_from_file(filename);
    //};

    //The stuff that I care about:
    PetscReal cos_factor() const { return cos_factor_; };
    PetscInt n_size() const { return n_size_; };
    PetscInt l_size() const { return l_size_; };
    PetscInt m_size() const { return m_size_; };

    //PetscErrorCode init_from_file();
    //PetscErrorCode save_parameters();

    std::string print() const;
private:
    void register_parameters();
    ez::ezOptionParser opt;
    int n_size_;
    int l_size_;
    int m_size_;
    double cos_factor_;
};

std::string AbsorberParameters::print() const
{
    std::ostringstream out;
    out << "absorber_n_size: " << n_size_ << std::endl;
    out << "absorber_l_size: " << l_size_ << std::endl;
    out << "absorber_m_size: " << m_size_ << std::endl;
    out << "absorber_cos_factor: " << cos_factor_ << std::endl;
    return out.str();
}

void AbsorberParameters::register_parameters()
{
    std::string prefix = "-absorber_";
    opt.overview = "Absorber Parameters";
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
            "40",
            0,
            1,
            0,
            "size of absorber in n",
            std::string(prefix).append("n_size\0").c_str()
           );
    opt.add(
            "10",
            0,
            1,
            0,
            "size of absorber in l",
            std::string(prefix).append("l_size\0").c_str()
           );
    opt.add(
            "0",
            0,
            1,
            0,
            "size of absorber in m",
            std::string(prefix).append("m_size\0").c_str()
           );
    opt.add(
            ".125",
            0,
            1,
            0,
            "the exponential on the cosine",
            std::string(prefix).append("cos_factor\0").c_str()
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
