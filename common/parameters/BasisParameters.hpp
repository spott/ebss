#pragma once
#include <unistd.h>
#include <common/common.hpp>
#include <common/parameters/Parameters.hpp>

#include <thread>

template <typename compute_type_, typename write_type_ = PetscReal>
class BasisParameters : public Parameters
{
  public:
    typedef compute_type_ compute_type;
    typedef write_type_ write_type;
    BasisParameters( int argc, const char** argv, MPI_Comm comm )
        : Parameters( comm ), procs_( 0 ), fs_( false ), l_only_( false ), n_only_( false ),
          charge_( 0 ), nmax_( 0 ), lmax_( 0 ), l_( 0 ), n_(0), points_( 0 ),
          rmax_( 0 ), rmin_( 0 ), bo_( false ), atom_( "hydrogen" )
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

        if ( opt.isSet( "-basis_config" ) ) {
            std::string fname;
            opt.get( "-basis_config" )->getString( fname );
            if ( !opt.importFile( fname.c_str(), '#' ) ) {
                std::cout << "file must exist!" << std::endl;
                throw std::exception();
            }
        }

        if ( opt.isSet( "-basis_l" ) ) {
            opt.get( "-basis_l" )->getInt( l_ );
            l_only_ = true;
        }
        if ( opt.isSet( "-basis_n" ) ) {
            opt.get( "-basis_n" )->getInt( n_ );
            n_only_ = true;
        }

        opt.get( "-basis_charge" )->getInt( charge_ );
        opt.get( "-basis_nmax" )->getInt( nmax_ );
        opt.get( "-basis_lmax" )->getInt( lmax_ );
        opt.get( "-basis_rmax" )->getDouble( rmax_ );
        opt.get( "-basis_rmin" )->getDouble( rmin_ );
        opt.get( "-basis_points" )->getInt( points_ );
        opt.get( "-basis_folder" )->getString( folder_ );
        opt.get( "-basis_atom" )->getString( atom_ );
        fs_ = opt.isSet( "-basis_fs" );
        bo_ = opt.isSet( "-basis_bound_only" );


        if ( opt.isSet( "-basis_procs" ) ) {
            // check for sign:
            int p = 0;
            opt.get( "-basis_procs" )->getInt( p );
            if ( p <= 0 )
                procs_ = std::thread::hardware_concurrency();
            else
                procs_ = static_cast<size_t>( p );
        } else
            procs_ = std::thread::hardware_concurrency();


        folder_ = common::absolute_path( folder_ );

        grid_ = std::vector<compute_type>( points_ );
        basis_prototype_ = std::vector<BasisID>();
    };

    BasisParameters( std::string filename, MPI_Comm comm )
        : Parameters( comm )
    {
        try { init_from_file( filename ); }
        catch ( std::exception e ) { throw e; }
    };

    // The stuff that I care about:
    PetscReal rmax() const
    {
        return rmax_;
    };
    PetscReal rmin() const
    {
        return rmin_;
    };
    PetscInt points() const
    {
        return points_;
    };
    PetscInt nmax() const
    {
        return nmax_;
    };
    PetscInt lmax() const
    {
        return lmax_;
    };
    PetscInt l() const
    {
        return l_;
    };
    PetscInt n() const
    {
        return n_;
    };
    PetscInt charge() const
    {
        return charge_;
    };
    bool fs() const
    {
        return fs_;
    };
    bool bound_only() const { return bo_; }
    bool l_only() const { return l_only_; }
    bool n_only() const { return n_only_; }
    size_t procs() const { return procs_; }
    std::string atom() const
    {
        return atom_;
    };

    // getting the folder:
    std::string basis_folder() const
    {
        return folder_;
    };
    std::string grid_filename() const
    {
        return std::string( folder_ ).append( "/grid.dat\0" );
    };
    std::string print() const;
    std::string basis_function_filename( BasisID a ) const;
    std::string l_block_filename( int l ) const;
    std::string basis_prototype_filename() const;
    void save_parameters();
    void init_from_file( std::string filename );

    std::vector<compute_type>& grid();
    std::vector<BasisID>& basis_prototype();

  private:
    std::vector<BasisID> basis_prototype_;
    std::vector<compute_type> grid_;

    void register_parameters();
    ez::ezOptionParser opt;
    std::string folder_;
    size_t procs_;
    bool fs_;
    bool l_only_;
    bool n_only_;
    int charge_;
    int nmax_;
    int lmax_;
    int l_;
    int n_;
    int points_;
    double rmax_;
    double rmin_;
    bool bo_;
    std::string atom_;
};

template <typename compute_type_, typename write_type_>
std::vector<compute_type_>&
BasisParameters<compute_type_, write_type_>::grid()
{
    return grid_;
}

template <typename compute_type_, typename write_type_>
std::vector<BasisID>&
BasisParameters<compute_type_, write_type_>::basis_prototype()
{
    return basis_prototype_;
}


template <typename compute_type_, typename write_type_>
void BasisParameters<compute_type_, write_type_>::init_from_file(
    std::string filename )
{
    register_parameters();
    opt.importFile( filename.c_str(), '#' );

    opt.get( "-basis_charge" )->getInt( charge_ );
    opt.get( "-basis_nmax" )->getInt( nmax_ );
    opt.get( "-basis_lmax" )->getInt( lmax_ );
    opt.get( "-basis_rmax" )->getDouble( rmax_ );
    opt.get( "-basis_rmin" )->getDouble( rmin_ );
    opt.get( "-basis_points" )->getInt( points_ );
    opt.get( "-basis_folder" )->getString( folder_ );
    opt.get( "-basis_atom" )->getString( atom_ );
    fs_ = opt.isSet( "-basis_fs" );
    bo_ = opt.isSet( "-basis_bound_only" );

    try
    {
        this->grid_ =
            common::vector_type_change<write_type_, compute_type_>(
                common::import_vector_binary<write_type_>(
                    this->grid_filename() ) );
    }
    catch ( std::exception e )
    {
        std::cout << "grid file not found" << std::endl;
        throw e;
    }

    try
    {
        this->basis_prototype_ = common::import_vector_binary<BasisID>(
            this->basis_prototype_filename() );
    }
    catch ( std::exception e )
    {
        std::cout << "prototype file not found" << std::endl;
        throw e;
    }
}

template <typename compute_type_, typename write_type_>
void BasisParameters<compute_type_, write_type_>::save_parameters()
{
    if ( this->rank() == 0 ) {
        common::export_vector_binary(
            grid_filename(),
            common::vector_type_change<compute_type_, write_type_>(
                this->grid_ ) );
        common::export_vector_binary( this->basis_prototype_filename(),
                                      this->basis_prototype_ );

        std::string f = std::string( folder_ ).append( "/Basis.config\0" );

        // write the file myself:
        std::ofstream file;
        file.open( f );
        file << "-basis_nmax " << nmax_ << std::endl;
        file << "-basis_lmax " << lmax_ << std::endl;
        file << "-basis_rmax " << rmax_ << std::endl;
        file << "-basis_rmin " << rmin_ << std::endl;
        file << "-basis_points " << points_ << std::endl;
        file << "-basis_folder " << folder_ << std::endl;
        file << "-basis_atom " << atom_ << std::endl;
        file << "-basis_charge " << charge_ << std::endl;
        if ( fs_ ) file << "-basis_fs" << std::endl;
        if ( bo_ ) file << "-basis_bound_only" << std::endl;
        file.close();
    }
}

template <typename compute_type_, typename write_type_>
std::string
BasisParameters<compute_type_, write_type_>::l_block_filename( int l )
    const
{
    std::ostringstream ss;
    ss << folder_;
    if ( fs_ )
        std::cerr << " not implemented yet " << std::endl;
    else
        ss << "/l_" << l << ".dat";
    return ss.str();
}

template <typename compute_type_, typename write_type_>
std::string
BasisParameters<compute_type_, write_type_>::basis_function_filename(
    BasisID a ) const
{
    std::ostringstream ss;
    ss << folder_;
    if ( fs_ )
        ss << "/n_" << a.n << "_l_" << a.l << "_j_" << a.j << ".dat";
    else
        ss << "/n_" << a.n << "_l_" << a.l << ".dat";
    return ss.str();
}

template <typename compute_type_, typename write_type_>
std::string
BasisParameters<compute_type_, write_type_>::basis_prototype_filename()
    const
{
    std::stringstream ss;
    ss << folder_;
    ss << "prototype.dat";
    return ss.str();
}

template <typename compute_type_, typename write_type_>
std::string BasisParameters<compute_type_, write_type_>::print() const
{
    std::ostringstream out;
    out << "basis_nmax: " << this->nmax_ << std::endl;
    out << "basis_lmax: " << lmax_ << std::endl;
    out << "basis_rmax: " << rmax_ << std::endl;
    out << "basis_rmin: " << rmin_ << std::endl;
    out << "basis_points: " << points_ << std::endl;
    out << "basis_folder: " << folder_ << std::endl;
    out << "basis_atom: " << atom_ << std::endl;
    out << "basis_procs: " << procs_ << std::endl;
    out << "basis_charge " << charge_ << std::endl;
    if ( fs_ ) out << "basis_fs" << std::endl;
    if ( bo_ ) out << "-basis_bound_only" << std::endl;
    return out.str();
}

template <typename compute_type_, typename write_type_>
void BasisParameters<compute_type_, write_type_>::register_parameters()
{
    std::string prefix = "-basis_";
    opt.overview = "Basis Parameters";
    opt.add( "", // Default.
             0,  // Required?
             0,  // Number of args expected.
             0,  // Delimiter if expecting multiple args.
             "Display usage instructions.", // Help description.
             "-h",                          // Flag token.
             "-help",                       // Flag token.
             "--help",                      // Flag token.
             "--usage"                      // Flag token.
             );
    opt.add( "0",
             0,
             1,
             0,
             "charge of the atom",
             std::string( prefix ).append( "charge\0" ).c_str() );
    opt.add( "",
             0,
             0,
             0,
             "finestructure included?",
             std::string( prefix ).append( "fs\0" ).c_str() );
    opt.add( "",
             0,
             0,
             0,
             "bound only?",
             std::string( prefix ).append( "bound_only\0" ).c_str() );
    opt.add( "",
             0,
             1,
             0,
             "number of processors",
             std::string( prefix ).append( "procs\0" ).c_str() );
    opt.add( "500",
             0,
             1,
             0,
             "Max n value",
             std::string( prefix ).append( "nmax\0" ).c_str() );
    opt.add( "50",
             0,
             1,
             0,
             "Max l value",
             std::string( prefix ).append( "lmax\0" ).c_str() );
    opt.add( "0",
             0,
             1,
             0,
             "single l value calculation (testing)",
             std::string( prefix ).append( "l\0" ).c_str() );
    opt.add( "1",
             0,
             1,
             0,
             "single n value calculation (testing)",
             std::string( prefix ).append( "n\0" ).c_str() );
    opt.add( "1000.",
             0,
             1,
             0,
             "Max r for grid",
             std::string( prefix ).append( "rmax\0" ).c_str() );
    opt.add( "0.000001",
             0,
             1,
             0,
             "Min r for grid",
             std::string( prefix ).append( "rmin\0" ).c_str() );
    opt.add( "10000",
             0,
             1,
             0,
             "Number of points on the grid",
             std::string( prefix ).append( "points\0" ).c_str() );
    opt.add( "hydrogen",
             1,
             1,
             0,
             "the atom to simulate",
             std::string( prefix ).append( "atom\0" ).c_str() );
    opt.add( "./",
             1,
             1,
             0,
             "Where the vectors are held",
             std::string( prefix ).append( "folder\0" ).c_str() );
    opt.add( "",
             0,
             1,
             0,
             "Config file to import",
             std::string( prefix ).append( "config\0" ).c_str() );
    opt.add(
        "", // Default.
        0,  // Required?
        0,  // Number of args expected.
        0,  // Delimiter if expecting multiple args.
        "Print all inputs and categories for debugging.", // Help
                                                          // description.
        "+d",
        "--debug" // Flag token.
        );
}
