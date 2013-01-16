#pragma once
#include<unistd.h>
#include<common/common.hpp>
#include<common/parameters/Parameters.hpp>
#include<common/parameters/BasisParameters.hpp>

template <typename write_type_ = double>
class HamiltonianParameters: public Parameters
{
public:
    HamiltonianParameters(int argc, const char** argv, MPI_Comm comm): Parameters(comm)
    {
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

        if (opt.isSet("-hamiltonian_config"))
        {
            std::string fname;
            opt.get("-hamiltonian_config")->getString(fname);
            if (! opt.importFile(fname.c_str(), '#'))
            {
                std::cout << "file must exist!" << std::endl;
                throw std::exception();
            }
        }

        opt.get("-hamiltonian_nmax")->getInt(nmax_);
        opt.get("-hamiltonian_lmax")->getInt(lmax_);
        opt.get("-hamiltonian_mmax")->getInt(mmax_);
        opt.get("-hamiltonian_folder")->getString(folder_);


        folder_ = common::absolute_path(folder_);

        opt.get("-hamiltonian_basis_config")->getString(basis_config_);
        basis_config_ = common::absolute_path(basis_config_);
        basis_ = new BasisParameters<write_type_, write_type_>(basis_config_, comm_);
        std::cout << basis_config_ << std::endl;

        //If we are bigger than the basis we are using, shrink to fit, and print an error message:
        if (this->nmax() > this->basis_->nmax() || this->lmax() > this->basis_->lmax() )
        {
            if (this->rank() == 0) 
            {
                std::cerr << "\e[31m======================ERROR========================" << std::endl;
                std::cerr << "=The basis you are attempting to use is smaller   =" << std::endl;
                std::cerr << "= than you wanted, change your values.            =" << std::endl;
                std::cerr << "======================ERROR========================\e[0m" << std::endl;
            }
            if (nmax_ > this->basis_->nmax())
                nmax_ = this->basis_->nmax();
            if (lmax_ > this->basis_->lmax())
                lmax_ = this->basis_->lmax();
        }

        //Initialize the prototype file:
        if (this->mmax() == 0 && this->nmax() == this->basis_->nmax() && this->lmax() == this->basis_->lmax())
        {
            this->prototype_.resize(0);
            for(auto i: *(this->basis_->basis_prototype()))
                this->prototype_.push_back(i);
            //std::copy(this->basis_->basis_prototype()->begin(), this->basis_->basis_prototype()->end(), this->prototype_);
            this->prototype_.shrink_to_fit();
        }
        else
        {
            BasisID temp;
            for(int l = 0; l <= this->lmax(); l++)
            {
                for (int j = 2 * l - 1 < 0 ? 2 * l + 1 : 2 * l - 1 ; j <= 2 * l + 1; j+=2)
                {
                    for(int m = -( l <= this->mmax() ? l : this->mmax() ); m <= ( l <= this->mmax() ? l : this->mmax() ); m++)
                    {
                        for (int n = 1; n <= this->nmax(); n++)
                        {
                            if (n > l)
                            {
                                temp.n = n;
                                temp.l = l;
                                temp.j = j;
                                auto el = std::find_if(this->basis_->basis_prototype()->begin(), 
                                        this->basis_->basis_prototype()->end(),
                                        [temp](BasisID a)->bool
                                        { return (a.n == temp.n && a.l == temp.l && a.j == temp.j);});
                                if ( el != this->basis_->basis_prototype()->end())
                                    temp.e = el->e;
                                else
                                {
                                    std::cout << "The energy wasn't in the basis prototype... I don't know what it is." << std::endl;
                                    throw(std::exception());
                                }
                                temp.m = m;
                                this->prototype_.push_back(temp);
                            }
                        }
                    }
                }
            }
            this->prototype_.shrink_to_fit();
        }
    };

    HamiltonianParameters(MPI_Comm comm, std::string fname): Parameters(comm)
    {
        init_from_file(fname);
    }

    int nmax() const { return nmax_; };
    int lmax() const { return lmax_; };
    int mmax() const { return mmax_; };

    std::string hamiltonian_folder() const { return folder_; };
    std::string dipole_matrix_filename() const;
    void write_dipole_matrix(Mat D);
    Mat read_dipole_matrix();
    std::string energy_eigenvalues_filename() const;
    void write_energy_eigenvalues();
    Vec read_energy_eigenvalues();
    std::string prototype_filename() const;
    BasisParameters<write_type_, write_type_>* basis_parameters() { return basis_; };

    std::vector<BasisID> prototype() const { return prototype_; };

    std::string print() const;

    void save_parameters();
    void init_from_file(std::string filename);

private:
    BasisParameters<write_type_, write_type_> *basis_;
    std::vector<BasisID> prototype_;

    void register_parameters();
    ez::ezOptionParser opt;
    int nmax_;
    int lmax_;
    int mmax_;
    std::string folder_;
    std::string basis_config_;
};

template< typename write_type_ >
void HamiltonianParameters<write_type_>::save_parameters()
{
    common::export_vector_binary(this->prototype_filename(), this->prototype_);
    //opt.exportFile(std::string(folder_).append("/Hamiltonian.config\0").c_str(), true);

    std::ofstream file;
    file.open(std::string(folder_).append("/Hamiltonian.config\0"));
    file << "-hamiltonian_nmax " << nmax_ << std::endl;
    file << "-hamiltonian_lmax " << lmax_ << std::endl;
    file << "-hamiltonian_mmax " << mmax_ << std::endl;
    file << "-hamiltonian_folder " << folder_ << std::endl;
    file << "-hamiltonian_basis_config " << basis_config_ << std::endl;
    file.close();
}

template< typename write_type_ >
void HamiltonianParameters<write_type_>::init_from_file(std::string filename)
{
    register_parameters();
    opt.importFile(filename.c_str(), '#');

    opt.get("-hamiltonian_nmax")->getInt(nmax_);
    opt.get("-hamiltonian_lmax")->getInt(lmax_);
    opt.get("-hamiltonian_mmax")->getInt(mmax_);
    opt.get("-hamiltonian_folder")->getString(folder_);
    opt.get("-hamiltonian_basis_config")->getString(basis_config_);
    basis_ = new BasisParameters<write_type_, write_type_>(basis_config_, comm_);

    this->prototype_ = common::import_vector_binary<BasisID>(this->prototype_filename());
}

template<typename write_type_ >
void HamiltonianParameters<write_type_>::write_dipole_matrix(Mat D)
{
    common::petsc_binary_write(this->dipole_matrix_filename(), D, this->comm_);
}

template<typename write_type_ >
Mat HamiltonianParameters<write_type_>::read_dipole_matrix()
{
    return common::petsc_binary_read<Mat>(this->dipole_matrix_filename(), this->comm_);
}

template<typename write_type_ >
void HamiltonianParameters<write_type_>::write_energy_eigenvalues()
{
    Vec ev;
    VecCreate(this->comm_, &ev);
    VecSetType(ev, VECSTANDARD);
    VecSetFromOptions(ev);
    VecSetSizes(ev, PETSC_DECIDE, this->prototype_.size());
    for(size_t i = 0; i < this->prototype_.size(); i++)
    {
        VecSetValue(ev, i, this->prototype_[i].e, INSERT_VALUES);
    }

    VecAssemblyBegin(ev);
    VecAssemblyEnd(ev);

    common::petsc_binary_write(this->energy_eigenvalues_filename(), ev, this->comm_);
}

template<typename write_type_ >
Vec HamiltonianParameters<write_type_>::read_energy_eigenvalues()
{
    return common::petsc_binary_read<Vec>(this->energy_eigenvalues_filename(), this->comm_);
}

template< typename write_type_ >
std::string HamiltonianParameters<write_type_>::dipole_matrix_filename() const
{
    std::stringstream ss;
    ss << folder_;
    ss << "/dipole_matrix.dat\0";
    return ss.str();
}


template< typename write_type_ >
std::string HamiltonianParameters<write_type_>::energy_eigenvalues_filename() const
{
    std::stringstream ss;
    ss << folder_;
    ss << "/energy_eigenvalues_vector.dat\0";
    return ss.str();
}


template< typename write_type_ >
std::string HamiltonianParameters<write_type_>::prototype_filename() const
{
    std::stringstream ss;
    ss << folder_;
    ss << "/vector_prototype.dat\0";
    return ss.str();
}

template<typename write_type_>
std::string HamiltonianParameters<write_type_>::print() const
{
    std::ostringstream out;
    out << "hamiltonian_nmax: " << nmax_ << std::endl;
    out << "hamiltonian_lmax: " << lmax_ << std::endl;
    out << "hamiltonian_mmax: " << mmax_ << std::endl;
    out << "hamiltonian_folder" << folder_ << std::endl;
    out << basis_->print();
    return out.str();
}

template< typename write_type_ >
void HamiltonianParameters<write_type_>::register_parameters()
{
    std::string prefix = "-hamiltonian_";
    opt.overview = "Hamiltonian Parameters";
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
            "500",
            0,
            1,
            0,
            "Max n value",
            std::string(prefix).append("nmax\0").c_str()
           );
    opt.add(
            "50",
            0,
            1,
            0,
            "Max l value",
            std::string(prefix).append("lmax\0").c_str()
           );
    opt.add(
            "0",
            0,
            1,
            0,
            "Max m value (abs)",
            std::string(prefix).append("mmax\0").c_str()
           );
    opt.add(
            "",
            1,
            1,
            0,
            "Basis config to import",
            std::string(prefix).append("basis_config\0").c_str()
           );
    opt.add(
            "./",
            0,
            1,
            0,
            "folder for hamiltonian",
            std::string(prefix).append("folder\0").c_str()
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
