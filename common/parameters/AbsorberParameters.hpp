#pragma once

#include<common/common.hpp>
#include<common/parameters/Parameters.hpp>

enum abs_type {
    COSINE, 
    CX_ROT,
    CX_SCALE
};

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

        std::string t;
        opt.get("-absorber_type")->getString(t);
        if (t == "cosine")
            type_ = COSINE;
        else if (t == "cx_rot")
            type_ = CX_ROT;
        else if (t == "cx_scale")
            type_ = CX_SCALE;
        else
        {
            std::cerr << "absorber type " << t << " not supported, defaulting to cosine" << std::endl;
            type_ = COSINE;
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
    std::string type() const {
        if (type_ == COSINE)
            return std::string("cosine");
        else if (type_ == CX_ROT)
            return std::string("cx_rot");
        else if (type_ == CX_SCALE)
            return std::string("cx_scale");
        else
            return std::string("no type!");
        }

    //create the absorber:
    void absorb(Vec *abs, HamiltonianParameters<PetscReal> *hparams);

    //PetscErrorCode init_from_file();
    void save_parameters();

    std::string print() const;
private:
    void register_parameters();
    ez::ezOptionParser opt;
    int n_size_;
    int l_size_;
    int m_size_;
    double cos_factor_;
    abs_type type_;
};

void AbsorberParameters::absorb(Vec *abs, HamiltonianParameters<PetscReal> *hparams)
{
    PetscInt start, end;
    VecGetOwnershipRange(*abs,&start,&end);
    VecSet(*abs, 1.);
    VecAssemblyBegin(*abs);
    VecAssemblyEnd(*abs);

    auto prototype = hparams->prototype();

    if (n_size_ == 0 && l_size_ == 0 && m_size_ == 0)
        return;
    
    if ( type_ == COSINE )
    {
        PetscReal val;
        if (rank() == 0)
            std::cerr << "cos_factor: " << cos_factor_ << std::endl;


        for (size_t i = start; i < end; i++)
        {
            val = 1.;
            if ((hparams->nmax() - prototype[i].n) < n_size())
                val *= std::pow(std::sin(
                            ((hparams->nmax() - prototype[i].n) * math::PI)/(2*n_size())
                            ),cos_factor());
            if ((hparams->lmax() - prototype[i].l) < l_size())
                val *= std::pow(std::sin(
                            ((hparams->lmax() - prototype[i].l) * math::PI)/(2*l_size())),
                        cos_factor());
            if ((hparams->mmax() - std::abs(prototype[i].m)) < m_size())
                val *= std::pow(std::sin(
                            ((hparams->mmax() - prototype[i].m) * math::PI)/(2*m_size())),
                        cos_factor());
            VecSetValue(*abs, i, val, INSERT_VALUES);
        }
    }
    else if ( type_ == CX_ROT )
    {
        PetscScalar val;
        std::cerr << "complex rotation is still being tested!  no l-absorber yet" << std::endl;
        for (size_t i = start; i < end; i++)
        {
            val = 1.;
            if ((hparams->nmax() - prototype[i].n) < n_size())
                val *= std::exp(
                            std::complex<double>(0,((n_size() - hparams->nmax() + prototype[i].n) * math::PI)/(2*n_size()))
                            );
            //if ((hparams->lmax() - prototype[i].l) < l_size())
                //val *= std::exp(
                            //std::complex<double>(0,((hparams->lmax() - prototype[i].l) * math::PI)/(2*l_size()))
                            //);
            //if ((hparams->mmax() - std::abs(prototype[i].m)) < m_size())
                //val *= std::exp(
                            //std::complex<double>(0,((hparams->mmax() - prototype[i].m) * math::PI)/(2*m_size()))
                            //);
            VecSetValue(*abs, i, val, INSERT_VALUES);
        }
    }

    VecAssemblyBegin(*abs);
    VecAssemblyEnd(*abs);

}

std::string AbsorberParameters::print() const
{
    std::ostringstream out;
    out << "absorber_n_size: " << n_size_ << std::endl;
    out << "absorber_l_size: " << l_size_ << std::endl;
    out << "absorber_m_size: " << m_size_ << std::endl;
    out << "absorber_cos_factor: " << cos_factor_ << std::endl;
    out << "absorber_type: " << type_ << std::endl;
    return out.str();
}

void AbsorberParameters::save_parameters()
{
    std::ofstream file;
    file.open(std::string("./Absorber.config"));
    file << "-absorber_n_size " << n_size_ << std::endl;
    file << "-absorber_l_size " << l_size_ << std::endl;
    file << "-absorber_m_size " << m_size_ << std::endl;
    file << "-absorber_cos_factor " << cos_factor_ << std::endl;
    if (type_ == COSINE)
        file << "-absorber_type cosine" << std::endl;
    if (type_ == CX_ROT)
        file << "-absorber_type cx_rot" << std::endl;
    if (type_ == CX_SCALE)
        file << "-absorber_type cx_scale" << std::endl;
    file.close();
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
            "cosine",
            0,
            1,
            0,
            "type (one of: cosine, cx_rot, or cx_scaled)",
            std::string(prefix).append("type\0").c_str()
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
