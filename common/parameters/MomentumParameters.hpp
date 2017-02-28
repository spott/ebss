#pragma once

#include <common/common.hpp>
#include <common/parameters/HamiltonianParameters.hpp>
#include <common/parameters/Parameters.hpp>
#include <common/types.hpp>


template <typename write_type_ = double>
class MomentumParameters : public Parameters
{
  public:
    MomentumParameters( int argc, const char** argv, MPI_Comm comm )
        : Parameters( comm )
    {
        register_parameters();

        opt.parse( argc, argv );

        if ( opt.isSet( "+d" ) ) {
            std::string pretty;
            opt.prettyPrint( pretty );
            std::cout << pretty;
        }
        if ( opt.isSet( "-h" ) ) {
            std::string usage;
            opt.getUsage( usage, 80, ez::ezOptionParser::ALIGN );
            std::cout << usage;
        }

        if ( opt.isSet( "-momentum_config" ) ) {
            std::string fname;
            opt.get( "-momentum_config" )->getString( fname );
            // init_from_file(fname);
            if ( !opt.importFile( fname.c_str(), '#' ) ) {
                std::cout << "file must exist!" << std::endl;
                throw std::exception();
            }
            opt.get( "-momentum_folder" )->getString( folder_ );
            folder_    = common::absolute_path( folder_ );
            kPrototype = common::import_vector_binary<kBasisID>(
                this->prototype_filename() );
        }

        opt.get( "-momentum_wf" )->getStrings( wf_filenames_ );

        opt.get( "-momentum_dtheta" )->getDouble( dtheta_ );
        opt.get( "-momentum_kmax" )->getDouble( kmax_ );
        opt.get( "-momentum_dk" )->getDouble( dk_ );
        opt.get( "-momentum_lmax" )->getInt( lmax_ );
        opt.get( "-momentum_out_filename_postfix" )->getString( ofnamepost_ );
        // opt.get("-momentum_folder")->getString(folder_);
        // folder_ = common::absolute_path(folder_);

        opt.get( "-momentum_nmax" )->getInt( nmax_ );
        opt.get( "-momentum_hamiltonian_config" )
            ->getString( hamiltonian_config_ );
        hamiltonian_config_ = common::absolute_path( hamiltonian_config_ );
        hamiltonian_        = new HamiltonianParameters<write_type_>(
            comm_, hamiltonian_config_ );
        // std::cout << hamiltonian_config_ << std::endl;

        std::cout << "wfs: " << std::endl;
        for ( auto& a : wf_filenames_ ) {
            std::cout << a << std::endl;
        }

        if ( nmax_ > hamiltonian_->nmax() ) nmax_ = hamiltonian_->nmax();
        if ( lmax_ > hamiltonian_->lmax() ) lmax_ = hamiltonian_->lmax();

        // full prototype:
        prototype_ = this->hamiltonian_->prototype();
        if ( kPrototype.size() == 0 ) kPrototype = genPrototype();
    }

    int    nmax() const { return nmax_; };
    int    lmax() const { return lmax_; };
    double kmax() const { return kmax_; };
    double dk() const { return dk_; };
    double dtheta() const { return dtheta_; };

    const std::vector<std::string>& wf_filenames() const
    {
        return wf_filenames_;
    }
    std::string output_filename_postfix() const { return ofnamepost_; };
    std::string momentum_folder() const { return folder_; };
    std::string matrix_filename() const
    {
        return std::string( momentum_folder() ).append( "/matrix.dat" );
    }
    std::string prototype_filename() const
    {
        return std::string( momentum_folder() ).append( "/kprototype.dat" );
    }
    std::vector<kBasisID>        genPrototype();
    const std::vector<kBasisID>& kprototype() const { return kPrototype; };
    // const std::vector< BasisID >& prototype() const { return prototype_; };
    const std::vector<BasisID>& prototype() const { return prototype_; };
    const HamiltonianParameters<write_type_>& hamiltonian() const
    {
        return *hamiltonian_;
    };
    void write_matrix( Mat ) const;
    void write_prototype() const;

    void save_parameters();
    void init_from_file( std::string filename );
    std::string print() const;

  private:
    // momentum quantum numbers:
    double kmax_;
    double dk_;
    int    lmax_;
    bool   fs = true;

    // pic:
    double dtheta_;

    // basis:
    int                                 nmax_;
    std::string                         hamiltonian_config_;
    HamiltonianParameters<write_type_>* hamiltonian_;

    // folder:
    std::string folder_;
    std::string ofnamepost_;

    void               register_parameters();
    ez::ezOptionParser opt;

    std::vector<kBasisID> kPrototype;
    // std::vector< BasisID > fullPrototype;
    std::vector<BasisID>     prototype_;
    std::vector<std::string> wf_filenames_;
};

template <typename write_type_>
void MomentumParameters<write_type_>::save_parameters()
{
    common::export_vector_binary( this->prototype_filename(),
                                  this->kPrototype );

    std::ofstream file;
    file.open( std::string( folder_ ).append( "/Momentum.config\0" ) );
    file << "-momentum_nmax " << nmax_ << std::endl;
    file << "-momentum_lmax " << lmax_ << std::endl;
    file << "-momentum_kmax " << kmax_ << std::endl;
    file << "-momentum_dk " << dk_ << std::endl;
    file << "-momentum_folder " << folder_ << std::endl;
    file << "-momentum_hamiltonian_config " << hamiltonian_config_ << std::endl;
    file.close();
}

template <typename write_type_>
void MomentumParameters<write_type_>::init_from_file( std::string filename )
{
    register_parameters();
    opt.importFile( filename.c_str(), '#' );

    opt.get( "-momentum_nmax" )->getInt( nmax_ );
    opt.get( "-momentum_lmax" )->getInt( lmax_ );
    opt.get( "-momentum_kmax" )->getDouble( kmax_ );
    opt.get( "-momentum_dk" )->getDouble( dk_ );
    opt.get( "-momentum_folder" )->getString( folder_ );
    opt.get( "-momentum_hamiltonian_config" )->getString( hamiltonian_config_ );
    hamiltonian_ =
        new HamiltonianParameters<write_type_>( comm_, hamiltonian_config_ );

    // this->prototype_ =
    // common::import_vector_binary<BasisID>(this->prototype_filename());
    kPrototype =
        common::import_vector_binary<kBasisID>( this->prototype_filename() );
    // kPrototype = genPrototype();
    std::cout << "inited from file" << std::endl;
}

template <typename write_type_>
std::string MomentumParameters<write_type_>::print() const
{
    std::ostringstream out;
    out << "momentum_nmax: " << nmax_ << std::endl;
    out << "momentum_lmax: " << lmax_ << std::endl;
    out << "momentum_kmax: " << kmax_ << std::endl;
    out << "momentum_dk: " << dk_ << std::endl;
    out << "momentum_dtheta: " << dtheta_ << std::endl;
    out << "momentum_folder" << folder_ << std::endl;
    for ( auto& a : wf_filenames_ ) out << "momentum_wf " << a << std::endl;
    out << hamiltonian_->print();
    return out.str();
}

template <typename write_type_>
void MomentumParameters<write_type_>::write_matrix( Mat D ) const
{
    common::petsc_binary_write( this->matrix_filename(), D, this->comm_ );
}

template <typename write_type_>
void MomentumParameters<write_type_>::write_prototype() const
{
    common::export_vector_binary( this->prototype_filename(),
                                  this->kPrototype );
}

// make the kPrototype:
template <typename write_type_>
std::vector<kBasisID> MomentumParameters<write_type_>::genPrototype()
{
    int size = ( static_cast<int>( kmax_ / dk_ ) + 1 ) * ( lmax_ + 1 );
    // std::cout << "size: " << size << std::endl;

    std::vector<kBasisID> kp;
    kp.reserve( size );
    double k = 0;

    for ( size_t i = 1; i <= kmax_ / dk_; i++ ) {
        k = i * dk_;
        for ( int l = 0; l <= lmax_; l++ ) kp.push_back( {k, l} );
    }

    kp.shrink_to_fit();

    // for (auto& a: kp)
    // std::cout << a.k << ", " << a.l << std::endl;

    return kp;
}

template <typename write_type_>
void MomentumParameters<write_type_>::register_parameters()
{
    std::string prefix = "-momentum_";
    opt.overview       = "Momentum Parameters";
    opt.add( "",  // Default.
             0,   // Required?
             0,   // Number of args expected.
             0,   // Delimiter if expecting multiple args.
             "Display usage instructions.",  // Help description.
             "-h",                           // Flag token.
             "-help",                        // Flag token.
             "--help",                       // Flag token.
             "--usage"                       // Flag token.
             );
    opt.add( "./wf_final.dat", 0, 1, 0, "filename of wf",
             std::string( prefix ).append( "wf\0" ).c_str() );
    opt.add( ".31415", 0, 1, 0, "the theta resolution",
             std::string( prefix ).append( "dtheta\0" ).c_str() );
    opt.add( "500", 0, 1, 0, "Max n value",
             std::string( prefix ).append( "nmax\0" ).c_str() );
    opt.add( "50", 0, 1, 0, "Max l value",
             std::string( prefix ).append( "lmax\0" ).c_str() );
    opt.add( ".1", 0, 1, 0, "Max dk value",
             std::string( prefix ).append( "dk\0" ).c_str() );
    opt.add( "2", 0, 1, 0, "Max k value",
             std::string( prefix ).append( "kmax\0" ).c_str() );
    opt.add( "", 1, 1, 0, "hamiltonian config to import",
             std::string( prefix ).append( "hamiltonian_config\0" ).c_str() );
    opt.add( "./", 0, 1, 0, "folder for hamiltonian",
             std::string( prefix ).append( "folder\0" ).c_str() );
    opt.add( "kspectrum.dat", 0, 1, 0, "output filename postfix",
             std::string( prefix ).append( "out_filename_postfix\0" ).c_str() );
    // opt.add(
    //"",
    // 0,
    // 0,
    // 0,
    //"the basis doesn't have a finestructure",
    // std::string(prefix).append("nofs\0").c_str()
    //);
    opt.add( "", 0, 1, 0, "Config file to import",
             std::string( prefix ).append( "config\0" ).c_str() );
    opt.add(
        "",  // Default.
        0,   // Required?
        0,   // Number of args expected.
        0,   // Delimiter if expecting multiple args.
        "Print all inputs and categories for debugging.",  // Help description.
        "+d",
        "--debug"  // Flag token.
        );
}
