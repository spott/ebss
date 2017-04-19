#pragma once

// ebss:
#include <common/common.hpp>
#include <common/math.hpp>
#include <common/parameters/Parameters.hpp>

// stl:
#include <numeric>
#include <sstream>
#include <string>

// petsc:
//#include<petsc.h>

class LaserParameters : public Parameters
{
  public:
    LaserParameters( int argc, const char** argv, MPI_Comm comm )
        : Parameters( comm )
    {
        register_parameters();

        opt.parse( argc, argv );

        get_parameters();
    };

    LaserParameters( MPI_Comm comm ) : Parameters( comm )
    {
        // when the parsing will be done elsewhere:
    }

    // call after parsing:
    void get_parameters();

    PetscReal   lambda() const;
    PetscReal   frequency() const;
    PetscReal   intensity() const;
    PetscReal   cep() const;
    PetscReal   cycles() const;
    PetscReal   dt() const;
    PetscReal   dt_after() const;
    PetscReal   t_after() const;
    std::string laser_filename() const;
    PetscReal   pulse_length() const;
    std::string shape() const { return laser_shape_; };
    PetscScalar envelope( PetscReal t, PetscReal t_start ) const;
    virtual PetscScalar efield( const PetscReal t,
                                const PetscReal phase ) const;
    virtual PetscScalar efield( const PetscReal t ) const;
    std::pair<std::vector<double>, std::vector<double>>
    read_efield( std::string filename ) const;

    virtual std::string print() const;

    void save_parameters() const;

  protected:
    ez::ezOptionParser opt;
    void               register_parameters();
    double             lambda_;
    double             intensity_;
    double             cep_;
    double             cycles_;
    double             dt_;
    double             dt_after_;
    double             t_after_;
    double             energy_;
    double             start_height_;
    std::string        laser_filename_;
    std::string        laser_shape_;
    bool               load_efield;
    int                laser_front_shape_;
    int                laser_back_shape_;
};

void LaserParameters::get_parameters()
{
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

    if ( opt.isSet( "-laser_config" ) ) {
        std::string fname;
        opt.get( "-laser_config" )->getString( fname );
        if ( !opt.importFile( fname.c_str(), '#' ) ) {
            std::cout << "laser config file must exist!" << std::endl;
            throw std::exception();
        }
    }
    if ( opt.isSet( "-laser_energy" ) )
        opt.get( "-laser_energy" )->getDouble( energy_ );
    else
        energy_ = -1.0;

    load_efield = opt.isSet( "-laser_filename" );

    opt.get( "-laser_lambda" )->getDouble( lambda_ );
    opt.get( "-laser_intensity" )->getDouble( intensity_ );
    opt.get( "-laser_cep" )->getDouble( cep_ );
    opt.get( "-laser_cycles" )->getDouble( cycles_ );
    opt.get( "-laser_dt" )->getDouble( dt_ );
    opt.get( "-laser_dt_after" )->getDouble( dt_after_ );
    opt.get( "-laser_t_after" )->getDouble( t_after_ );
    opt.get( "-laser_filename" )->getString( laser_filename_ );
    opt.get( "-laser_front_shape" )->getInt( laser_front_shape_ );
    opt.get( "-laser_back_shape" )->getInt( laser_back_shape_ );
    opt.get( "-laser_envelope" )->getString( laser_shape_ );
    opt.get( "-laser_height" )->getDouble( start_height_ );
}

std::string LaserParameters::print() const
{
    std::ostringstream out;
    out << std::scientific;
    out.precision( 15 );
    out << "laser_lambda: " << lambda_ << std::endl;
    out << "laser_intensity: " << intensity_ << std::endl;
    out.unsetf( std::ios::floatfield );
    out.precision( 5 );
    out << "laser_cep: " << cep_ << std::endl;
    out << "laser_cycles: " << cycles_ << std::endl;
    out << "laser_dt: " << dt_ << std::endl;
    out << "laser_dt_after: " << dt_after_ << std::endl;
    out << "laser_t_after: " << t_after_ << std::endl;
    out << "laser_filename: " << laser_filename_ << std::endl;
    out << "laser_front_shape " << laser_front_shape_ << std::endl;
    out << "laser_back_shape " << laser_back_shape_ << std::endl;
    out << "laser_envelope " << laser_shape_ << std::endl;
    out << "laser_height " << start_height_ << std::endl;

    return out.str();
}
void LaserParameters::save_parameters() const
{
    std::ofstream file;
    file.open( std::string( "./Laser.config" ) );
    file << std::scientific;
    file.precision( 15 );
    file << "-laser_lambda " << lambda_ << std::endl;
    file << "-laser_intensity " << intensity_ << std::endl;
    file.unsetf( std::ios::floatfield );
    file.precision( 5 );
    file << "-laser_cep " << cep_ << std::endl;
    file << "-laser_cycles " << cycles_ << std::endl;
    file << "-laser_dt " << dt_ << std::endl;
    file << "-laser_dt_after " << dt_after_ << std::endl;
    file << "-laser_t_after " << t_after_ << std::endl;
    file << "-laser_filename " << laser_filename_ << std::endl;
    file << "-laser_front_shape " << laser_front_shape_ << std::endl;
    file << "-laser_back_shape " << laser_back_shape_ << std::endl;
    file << "-laser_envelope " << laser_shape_ << std::endl;
    file << "-laser_height " << start_height_ << std::endl;
    file.close();
}

PetscReal LaserParameters::lambda() const
{
    return ( lambda_ / 5.29177206e-2 );
}
PetscReal LaserParameters::frequency() const
{
    if ( energy_ < 0 )
        return ( math::C * 2 * math::PI ) / ( this->lambda() );
    else
        return energy_;
}
PetscReal LaserParameters::intensity() const
{
    return intensity_ / ( 3.5094452e16 );
}
//{ return intensity_/3.5101e+16; }
PetscReal LaserParameters::cep() const { return cep_; }
PetscReal LaserParameters::cycles() const { return cycles_; }
PetscReal LaserParameters::dt() const { return dt_; }
PetscReal LaserParameters::dt_after() const { return dt_after_; }
PetscReal LaserParameters::t_after() const { return t_after_; }
PetscReal LaserParameters::pulse_length() const
{
    if ( this->shape() == "sin_squared" )
        return ( ( math::PI * 2 * cycles() ) / frequency() ) + t_after();
    else if ( this->shape() == "gaussian" ) {
        PetscReal fwhm_time = ( math::PI * 2 * cycles() ) / frequency();
        PetscReal mean      = fwhm_time *
                         std::sqrt( std::log( 1. / this->start_height_ ) ) /
                         ( 2. * std::sqrt( std::log( 2. ) ) );
        return mean * 2.;
    }
}
std::string LaserParameters::laser_filename() const
{
    return std::string( laser_filename_ );
}

PetscScalar LaserParameters::envelope( PetscReal t, PetscReal t_start ) const
{
    if ( this->shape() == "sin_squared" ) {
        if ( t < t_start || t > t_start + ( pulse_length() - t_after() ) )
            return 0;

        PetscReal efield = 0.0;
        if ( t < t_start + pulse_length() / 2. )
            efield = std::pow( std::sin( this->frequency() * ( t - t_start ) /
                                         ( this->cycles() * 2 ) ),
                               laser_front_shape_ );
        if ( t >= t_start + pulse_length() / 2. )
            efield = std::pow( std::sin( this->frequency() * ( t - t_start ) /
                                         ( this->cycles() * 2 ) ),
                               laser_back_shape_ );
        return efield;
    } else if ( this->shape() == "gaussian" ) {
        PetscReal fwhm_time = ( math::PI * 2 * cycles() ) / frequency();
        PetscReal mean      = fwhm_time *
                             std::sqrt( std::log( 1. / this->start_height_ ) ) /
                             ( 2. * std::sqrt( std::log( 2. ) ) ) +
                         t_start;
        PetscReal std_deviation = fwhm_time / std::sqrt( 8. * std::log( 2. ) );

        return std::exp( -( t - mean ) * ( t - mean ) /
                         ( 2. * std_deviation * std_deviation ) );
    }
}

std::pair<std::vector<double>, std::vector<double>>
    LaserParameters::read_efield( std::string /*filename*/ ) const
{
    // ignore reading it in this time:
    std::vector<double> t( 100 );
    std::iota( t.begin(), t.end(), 0 );
    for ( auto& i : t ) i *= .01;

    std::vector<double> e( 100 );
    PetscReal           efield = std::sqrt( this->intensity() );
    for ( size_t i = 0; i < 100; ++i ) {
        e[i] = efield * std::pow( std::sin( this->frequency() * t[i] /
                                            ( this->cycles() * 2 ) ),
                                  2 ) *
               std::sin( this->frequency() * t[i] );
    }

    return std::make_pair( t, e );
}

PetscScalar LaserParameters::efield( PetscReal t, PetscReal phase ) const
{
    if ( t > pulse_length() || t < 0 ) return 0.0;
    // if ( t * this->frequency() / ( this->cycles() * 2 ) > math::PI || t
    // < 0 )
    PetscReal efield = std::sqrt( this->intensity() );
    return std::complex<double>(
        efield * envelope( t, 0.0 ) *
        // std::pow( std::sin( this->frequency() * t / ( this->cycles() * 2
        // ) ),
        // 2 ) *
        std::sin( this->frequency() * ( t - pulse_length()/2. ) + phase ) );
}


PetscScalar LaserParameters::efield( PetscReal t ) const
{
    if ( load_efield == false ) {
        return this->efield( t, this->cep() );
    } else {
        static math::interpolate efield = [=]() {
            auto a = read_efield( this->laser_filename_ );
            return math::interpolate( std::get<0>( a ), std::get<1>( a ) );
        }();
        return efield( t );
    }
}


void LaserParameters::register_parameters()
{
    std::string prefix = "-laser_";
    opt.overview       = "Laser Parameters";
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
    opt.add( "2", 0, 1, 0, "front shape",
             std::string( prefix ).append( "front_shape\0" ).c_str() );
    opt.add( "2", 0, 1, 0, "back shape",
             std::string( prefix ).append( "back_shape\0" ).c_str() );
    opt.add( "800", 0, 1, 0, "wavelength of the laser in nm",
             std::string( prefix ).append( "lambda\0" ).c_str() );
    opt.add( "", 0, 1, 0,
             "energy of the laser.  If this and the wavelength are specified, "
             "this will take precedence",
             std::string( prefix ).append( "energy\0" ).c_str() );
    opt.add( "sin_squared", 0, 1, 0, "Pulse envelope shape",
             std::string( prefix ).append( "envelope\0" ).c_str() );
    opt.add( "1e-16", 0, 1, 0,
             "for a gaussian pulse, the initial envelope value",
             std::string( prefix ).append( "height\0" ).c_str() );
    opt.add( "10e12", 0, 1, 0, "intensity of laser in W/cm^2",
             std::string( prefix ).append( "intensity\0" ).c_str() );
    opt.add( "0", 0, 1, 0, "Carrier Envelope Phase",
             std::string( prefix ).append( "cep\0" ).c_str() );
    opt.add( "10", 0, 1, 0, "number of cycles of the carrier frequency",
             std::string( prefix ).append( "cycles\0" ).c_str() );
    opt.add( "0.01", 0, 1, 0, "timestep durring the pulse (in A.U.)",
             std::string( prefix ).append( "dt\0" ).c_str() );
    opt.add( "0.01", 0, 1, 0, "timestep after the pulse (in A.U.)",
             std::string( prefix ).append( "dt_after\0" ).c_str() );
    opt.add( "0", 0, 1, 0, "time after the pulse (in A.U.)",
             std::string( prefix ).append( "t_after\0" ).c_str() );
    opt.add( "./efield.dat", 0, 1, 0,
             "filename to put the laser or pull the laser",
             std::string( prefix ).append( "filename\0" ).c_str() );
    opt.add( "", 0, 1, 0, "Config file to import",
             std::string( prefix ).append( "config\0" ).c_str() );
    opt.add( "",  // Default.
             0,   // Required?
             0,   // Number of args expected.
             0,   // Delimiter if expecting multiple args.
             "Print all inputs and categories for debugging.",  // Help
             // description.
             "+d",
             "--debug"  // Flag token.
             );
}
