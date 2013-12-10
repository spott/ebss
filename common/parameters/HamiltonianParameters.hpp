#pragma once
#include<common/common.hpp>
#include<common/output.hpp>
#include<common/parameters/Parameters.hpp>
#include<common/parameters/BasisParameters.hpp>
#include<algorithm>

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
            std::cout << "reading from file: ";
            std::string fname;
            opt.get("-hamiltonian_config")->getString(fname);
            std::cout << fname << std::endl;
            if (! opt.importFile(fname.c_str(), '#'))
            {
                std::cerr << "file must exist!" << std::endl;
                throw std::exception();
            }
        }

        //std::string ctemp;
        //if ( opt.isSet("-hamiltonian_nofs") )
        //{
            //fs = false;
            //std::cout << "no finestructure! "<< std::endl;
        //}
        //else
            //fs = true;
        //opt.get("-hamiltonian_components")->getString(ctemp);
        //if ( ctemp.find('x') != std::string::npos )
            //x = true;
        //if ( ctemp.find('y') != std::string::npos )
            //y = true;
        //if ( ctemp.find('z') != std::string::npos )
            //z = true;

        int m;
        opt.get("-hamiltonian_mem_per_proc")->getInt(m);

        if (m < 0)
            throw std::out_of_range("mem_per_proc is less than zero");
        else 
            mem_ = static_cast<size_t>(m)* 1024;


        opt.get("-hamiltonian_nmax")->getInt(nmax_);
        opt.get("-hamiltonian_lmax")->getInt(lmax_);
        opt.get("-hamiltonian_mmax")->getInt(mmax_);
        opt.get("-hamiltonian_folder")->getString(folder_);

        analytic_ = opt.isSet("-hamiltonian_analytical");

        folder_ = common::absolute_path(folder_);

        if (!analytic_)
        {
            opt.get("-hamiltonian_basis_config")->getString(basis_config_);
            basis_config_ = common::absolute_path(basis_config_);
            try {
            basis_ = new BasisParameters<write_type_, write_type_>(basis_config_, comm_);
            } catch (std::exception e) {
                std::cerr << "couldn't create basisparams" << std::endl;
            }

        //std::cout << basis_config_ << std::endl;
            fs_ = basis_->fs();
        }
        else
            basis_ = nullptr;

        //If we are bigger than the basis we are using, shrink to fit, and print an error message:
        //create the prototype file:

        this->prototype_.resize(0);
        if (!analytic_)
        {
            if (opt.isSet("-hamiltonian_config"))
                this->prototype_ = common::import_vector_binary<BasisID>(this->prototype_filename());
            else {
                this->prototype_.reserve(this->basis_->basis_prototype()->size());

                if (this->nmax() > this->basis_->nmax() || this->lmax() > this->basis_->lmax() )
                {
                    if (this->rank() == 0) 
                    {
                        std::cerr << output::red;
                        std::cerr << "======================ERROR========================" << std::endl;
                        std::cerr << "=The basis you are attempting to use is smaller   =" << std::endl;
                        std::cerr << "= than requested, your requested values have been =" << std::endl;
                        std::cerr << "= changed to fit                                  =" << std::endl;
                        std::cerr << "======================ERROR========================" << std::endl;
                        std::cerr << output::reset;
                    }
                    if (nmax_ > this->basis_->nmax())
                        nmax_ = this->basis_->nmax();
                    if (lmax_ > this->basis_->lmax())
                        lmax_ = this->basis_->lmax();
                }
                //make the prototype... be sure to order by l's first, then j's then n's.

                auto start = this->basis_->basis_prototype()->begin();
                auto end = this->basis_->basis_prototype()->end();
                if (fs_)
                    for(int l = 0; l <= lmax_; ++l )
                    {
                        for (int j = 2 * l - 1 ; j <= 2*l+1; j+=2)
                        {
                            if (j < 0)
                                continue;
                            for (int n = 1; n <= nmax_; ++n)
                            {
                                if (l <= n-1)
                                {
                                    auto loc = std::find_if(
                                            start, end, [n, l, j](const BasisID& a ) {
                                            return a.n == n && a.l == l && a.j == j;
                                            } );
                                    if (loc == end)
                                        throw std::out_of_range("didn't find the value we were looking for when making the prototype");
                                    else
                                        prototype_.push_back(*loc);
                                }
                                else continue;
                            }
                        }
                    }
                else
                    for(int l = 0; l <= lmax_; ++l )
                    {
                        for (int n = 1; n <= nmax_; ++n)
                        {
                            if (l <= n-1)
                            {
                                auto loc = std::find_if(
                                    start, end, [n, l](const BasisID& a ) {
                                        return a.n == n && a.l == l;
                                    } );
                                if (loc == end)
                                    throw std::out_of_range("didn't find the value we were looking for when making the prototype");
                                else
                                    prototype_.push_back(*loc);
                            }
                            else continue;
                        }
                    }
                //for(const auto i: *(this->basis_->basis_prototype()))
                //{
                    //if (i.n < this->nmax() && i.l < this->lmax())
                        //this->prototype_.push_back(i);
                //}
                this->prototype_.shrink_to_fit();
            }
        }
        else
        {
            //make the prototype:
            //go though l's first, then n's to reduce communication:
            for(int l = 0; l <= lmax_; ++l )
            {
                for (int n = 1; n <= nmax_; ++n)
                {
                    if (l <= n-1)
                        prototype_.push_back( {n, l, 0, 0, std::complex<double>( -1./(2. * n * n), 0) } );
                    else
                        continue;
                }
            }
        }
    }


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
    BasisParameters<write_type_, write_type_>* basis_parameters() const { return basis_; };

    std::vector<BasisID>& prototype() { return prototype_; };

    std::string print() const;
    bool fs() const { return fs_; }
    bool analytic() const { return analytic_; }
    size_t  mem_per_proc() const { return mem_; }

    void save_parameters();
    void init_from_file(std::string filename);

private:
    BasisParameters<write_type_, write_type_> *basis_;
    std::vector<BasisID> prototype_;

    bool fs_;
    void register_parameters();
    ez::ezOptionParser opt;
    int nmax_;
    int lmax_;
    int mmax_;
    size_t mem_;
    std::string folder_;
    std::string basis_config_;
    bool analytic_;
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
    if (!analytic_)
        file << "-hamiltonian_basis_config " << basis_config_ << std::endl;
    if (analytic_)
        file << "-hamiltonian_analytical" << std::endl;
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
    analytic_ = opt.isSet("-hamiltonian_analytical");
    if (!analytic_)
    {
        opt.get("-hamiltonian_basis_config")->getString(basis_config_);
        basis_ = new BasisParameters<write_type_, write_type_>(basis_config_, comm_);
    }
    else
        basis_ = nullptr;

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
    out << "hamiltonian_folder " << folder_ << std::endl;
	out << "hamiltonian_mem_per_proc " << mem_ << std::endl;
    if (!analytic_)
        out << basis_->print();
    else
        out << "hamiltonian_analytical" << std::endl;
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
            "",
            0,
            0,
            0,
            "calculate the analytical solution",
            std::string(prefix).append("analytical\0").c_str()
           );
    opt.add(
            "193274000",
            0,
            1,
            0,
            "mememory available per processor (in kb)",
            std::string(prefix).append("mem_per_proc\0").c_str()
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
            0,
            0,
            "the basis doesn't have a finestructure",
            std::string(prefix).append("nofs\0").c_str()
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
