#pragma once

// ebss:
#include <common/common.hpp>
#include <common/math.hpp>
#include <common/parameters/Parameters.hpp>

// stl:
#include <sstream>
#include <string>

// petsc:
//#include<petsc.h>

class StateParameters : public Parameters
{
  public:
    StateParameters( int argc, const char** argv, MPI_Comm comm )
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

        if ( opt.isSet( "-state_config" ) ) {
            std::string fname;
            opt.get( "-state_config" )->getString( fname );
            if ( !opt.importFile( fname.c_str(), '#' ) ) {
                std::cout << "state config file must exist!" << std::endl;
                throw std::exception();
            }
        }

        if ( opt.isSet( "-state_load" ) ) {
            std::string fn;
            opt.get( "-state_load" )->getString( fn );
            empty_states_ = common::import_vector_ascii<BasisID>( fn );
        } else
            empty_states_ = std::vector<BasisID>();

        if ( opt.isSet( "-state_no_bound" ) )
            nobound = true;
        else
            nobound = false;
        if ( opt.isSet( "-state_no_continuum" ) )
            nocontinuum = true;
        else
            nocontinuum = false;

        if ( opt.isSet( "-state_no_continuum_transitions" ) )
            notransitions = true;
        else
            notransitions = false;

        if ( opt.isSet( "-state_no_bound_continuum_transitions" ) )
            noboundcontinuum = true;
        else
            noboundcontinuum = false;

        opt.get( "-state_m" )->getInt( m_ );

        std::vector<std::vector<int>> init;
        opt.get( "-state_init" )->getMultiInts( init );
        init_.n = init[0][0];
        init_.l = init[0][1];
        init_.j = init[0][2];

        if ( opt.isSet( "-state_initial_wavefunction" ) )
            opt.get( "-state_initial_wavefunction" )
                ->getString( wavefunction_fname_ );


        opt.get( "-state_filename" )->getString( filename_ );
        filename_ = common::absolute_path( filename_ );

        std::vector<std::vector<int>> added_states;
        std::vector<std::vector<int>> removed_states;
        opt.get( "-state_add" )->getMultiInts( added_states );
        opt.get( "-state_rem" )->getMultiInts( removed_states );

        for ( size_t i = 0; i < removed_states.size(); i++ )
            empty_states_.push_back(
                {removed_states[i][0], removed_states[i][1], 0,
                 removed_states[i][2], std::complex<double>( 0 )} );
        for ( size_t i = 0; i < added_states.size(); i++ )
            add_states.push_back( {added_states[i][0], added_states[i][1], 0,
                                   added_states[i][2],
                                   std::complex<double>( 0 )} );
    };

    std::string print() const;
    void        save_parameters() const;

    std::vector<int> empty_states_index( const std::vector<BasisID> prototype );
    void disallowed_transitions( Mat&                        D,
                                 const std::vector<BasisID>& prototype );
    std::vector<BasisID> empty_states( const std::vector<BasisID> prototype );
    void initial_vector( Vec* v, const std::vector<BasisID> prototype );

  private:
    bool                 nobound;
    bool                 nocontinuum;
    bool                 notransitions, noboundcontinuum;
    ez::ezOptionParser   opt;
    void                 register_parameters();
    BasisID              init_;
    int                  m_;
    std::vector<BasisID> empty_states_;
    std::vector<BasisID> add_states;
    std::string          filename_;
    std::string          wavefunction_fname_;
    std::vector<std::array<int, 2>> matrix_elements_;
};

std::vector<BasisID>
StateParameters::empty_states( const std::vector<BasisID> prototype )
{
    if ( !nobound && !nocontinuum && m_ == 0 ) return empty_states_;

    for ( auto p : prototype ) {
        if ( p.e.real() < 0 && nobound &&
             !( p.n == init_.n && p.l == init_.l && p.j == init_.j ) ) {
            empty_states_.push_back( p );

        } else if ( p.e.real() > 0 && nocontinuum ) {
            empty_states_.push_back( p );
        } else if ( p.l < std::abs( p.m ) ) {
            std::cout << "removing state: " << std::endl;
            empty_states_.push_back( p );
        }
    }

    // auto add_states_it = add_states.begin();

    for ( auto a : add_states ) {
        auto e = empty_states_.begin();
        for ( ; e < empty_states_.end(); e++ ) {
            if ( ( ( *e ).n == a.n ) && ( ( *e ).l == a.l ) &&
                 ( ( *e ).j == a.j ) )
                empty_states_.erase( e );
        }
    }
    nobound     = false;
    nocontinuum = false;
    return empty_states_;
}

void StateParameters::disallowed_transitions(
    Mat& D, const std::vector<BasisID>& prototype )
{
    if ( !notransitions && !noboundcontinuum ) return;
    // std::vector<int> from_states;
    // std::vector<int> to_states;

    auto zero  = std::complex<double>( 0.0, 0.0 );
    int  total = 0;

    for ( auto i = prototype.cbegin(); i < prototype.cend(); ++i ) {
        for ( auto j = prototype.cbegin(); j < prototype.cend(); ++j ) {
            if ( notransitions )
                if ( j->e.real() > 0.0 && i->e.real() > 0.0 &&
                     ( std::abs( j->l - i->l ) == 1 ) ) {
                    std::vector<int>                  column;
                    std::vector<std::complex<double>> zeros;
                    while ( j->e.real() > 0.0 &&
                            ( std::abs( j->l - i->l ) == 1 ) &&
                            j < prototype.cend() ) {
                        column.push_back( j - prototype.cbegin() );
                        zeros.push_back( zero );
                        j++;
                    }
                    if ( rank() == 0 )
                        std::cout << "(" << *i << "->" << *( j - 1 ) << ") "
                                  << column.size() << std::endl;
                    int a = ( i - prototype.cbegin() );
                    MatSetValues( D, 1, &( a ), column.size(), column.data(),
                                  zeros.data(), INSERT_VALUES );
                }
            if ( noboundcontinuum )
                if ( i->e.real() < 0.0 && j->e.real() > 0.0 &&
                     ( std::abs( j->l - i->l ) == 1 ) ) {
                    std::vector<int>                  column;
                    std::vector<std::complex<double>> zeros;
                    while ( j->e.real() > 0.0 &&
                            ( std::abs( j->l - i->l ) == 1 ) &&
                            j < prototype.cend() ) {
                        column.push_back( j - prototype.cbegin() );
                        zeros.push_back( zero );
                        j++;
                    }
                    if ( rank() == 0 )
                        std::cout << "(" << *i << "->" << *( j - 1 ) << ") "
                                  << column.size() << std::endl;
                    int a = ( i - prototype.cbegin() );
                    MatSetValues( D, 1, &( a ), column.size(), column.data(),
                                  zeros.data(), INSERT_VALUES );
                }
        }
        MatAssemblyBegin( D, MAT_FLUSH_ASSEMBLY );
        MatAssemblyEnd( D, MAT_FLUSH_ASSEMBLY );
    }

    // return std::array<std::vector<int>, 2>{{from_states, to_states}};
}

std::vector<int>
StateParameters::empty_states_index( const std::vector<BasisID> prototype )
{
    empty_states( prototype );

    std::vector<int> state_index;
    for ( auto a : empty_states_ ) {
        int  i;
        auto it =
            std::find_if( prototype.begin(), prototype.end(), [a]( BasisID b ) {
                return ( a.n == b.n && a.l == b.l && a.j == b.j );
            } );
        if ( it == prototype.end() )
            std::cerr << a << " wasn't found in prototype" << std::endl;
        else {
            i = it - prototype.begin();
            state_index.push_back( i );
        }
    }

    return state_index;
}

void StateParameters::initial_vector( Vec*                       v,
                                      const std::vector<BasisID> prototype )
{
    if ( !wavefunction_fname_.empty() ) {
        if ( !common::file_exists( wavefunction_fname_ ) ) {
            std::cerr << wavefunction_fname_ << " doesn't exist" << std::endl;
            throw( std::exception() );
        }
        PetscViewer view;
        PetscViewerBinaryOpen( this->comm_, wavefunction_fname_.c_str(),
                               FILE_MODE_READ, &view );
        VecLoad( *v, view );
        PetscViewerDestroy( &view );
    } else {
        for ( size_t i = 0; i < prototype.size(); i++ ) {
            if ( prototype[i].n == init_.n && prototype[i].l == init_.l &&
                 prototype[i].j == init_.j ) {
                VecSetValue( *v, i, 1., INSERT_VALUES );
                break;
            }
        }
        VecAssemblyBegin( *v );
        VecAssemblyEnd( *v );
    }
}

std::string StateParameters::print() const
{
    std::ostringstream out;

    out << "state_no_bound " << nobound << std::endl;
    out << "state_no_continuum " << nocontinuum << std::endl;
    out << "state_filename " << filename_ << std::endl;
    out << "state_init " << init_ << std::endl;
    out << "state_m " << m_ << std::endl;

    return out.str();
}

// TODO fill this out
void StateParameters::save_parameters() const
{
    common::export_vector_ascii( filename_, empty_states_ );

    std::ofstream file;
    file.open( std::string( "./State.config" ) );
    if ( wavefunction_fname_.empty() )
        file << "-state_init " << init_.n << "," << init_.l << "," << init_.j
             << std::endl;
    else
        file << "-state_initial_wavefunction " << wavefunction_fname_
             << std::endl;
    for ( auto a : empty_states_ ) {
        file << "-state_rem " << a.n << "," << a.l << "," << a.j << std::endl;
    }
    file << "-state_m " << m_ << std::endl;
    file.close();
}

void StateParameters::register_parameters()
{
    std::string prefix = "-state_";
    opt.overview       = "State Parameters";

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

    opt.add( "", 0, 3, ',',
             "add a specific state (n,l,j triplet), or set of states (if "
             "removed otherwise)",
             std::string( prefix ).append( "add\0" ).c_str() );

    opt.add( "", 0, 3, ',',
             "remove a specific state, or set of states (n,l,j pair)",
             std::string( prefix ).append( "rem\0" ).c_str() );

    opt.add(
        "", 0, 0, 0, "remove continuum transitions",
        std::string( prefix ).append( "no_continuum_transitions\0" ).c_str() );

    opt.add( "", 0, 0, 0, "remove bound-continuum transitions",
             std::string( prefix )
                 .append( "no_bound_continuum_transitions\0" )
                 .c_str() );

    opt.add( "./empty_states.dat", 0, 0, 0, "load from file",
             std::string( prefix ).append( "load\0" ).c_str() );

    opt.add( "1,0,1", 0, 3, ',', "initial state: (n,l,j pair)",
             std::string( prefix ).append( "init\0" ).c_str() );

    opt.add( "", 0, 1, ',', "m state",
             std::string( prefix ).append( "m\0" ).c_str() );


    opt.add( "", 0, 1, 0,
             "initial wavefunction filename (PetscVec binary file)",
             std::string( prefix ).append( "initial_wavefunction\0" ).c_str() );

    opt.add( "", 0, 0, 0, "remove the bound states (toggle)",
             std::string( prefix ).append( "no_bound\0" ).c_str() );

    opt.add( "", 0, 0, 0, "remove the continuum states (toggle)",
             std::string( prefix ).append( "no_continuum\0" ).c_str() );

    opt.add( "./empty_states.dat", 0, 1, 0, "filename for states file",
             std::string( prefix ).append( "filename\0" ).c_str() );

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
