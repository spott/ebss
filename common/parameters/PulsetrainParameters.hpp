#pragma once

// ebss:
#include <common/common.hpp>
#include <common/math.hpp>
#include <common/parameters/LaserParameters.hpp>

// stl:
#include <sstream>
#include <string>

// petsc:
//#include<petsc.h>

class PulsetrainParameters : public LaserParameters
{
  public:
    PulsetrainParameters( int argc, const char** argv, MPI_Comm comm )
        : LaserParameters( argc, argv, comm )
    {
        register_parameters();

        opt.parse( argc, argv );

        get_parameters();
        s_au = spacing();

        if (test_)
        {
            auto a = test();
            common::export_vector_binary("./time_test.dat", a[0]);
            common::export_vector_binary("./efield_test.dat", a[1]);
        }
    }
    ;


    // call after parsing
    void get_parameters();

    bool in_pulse( PetscReal t ) const;
    std::vector<double> spacing() const;
    std::vector<double> phase() const;
    virtual PetscScalar efield( PetscReal t ) const;
    PetscReal max_time() const;
    PetscReal next_pulse_start( PetscReal t ) const;

    void save_parameters() const;
    std::string print() const;

    std::array<std::vector<double>,2> test() const;

  protected:
    bool test_;
    bool single_pulse;
    void register_parameters();
    std::vector<double> spacing_;
    std::vector<double> phase_;
    std::vector<double> s_au; // spacing in au
};

std::array<std::vector<double>, 2> PulsetrainParameters::test() const
{
    std::array<std::vector<double>, 2> out{
        {std::vector<double>( max_time() / this->dt() ),
         std::vector<double>( max_time() / this->dt() )}};


    out[0].clear();
    out[1].clear();
    //out[0].reserve( max_time() / this->dt() );
    //out[1].reserve( max_time() / this->dt() );

    double t = 0.0;

    while ( t < max_time() )
    {
        //if ( in_pulse(t) )
        //{
        out[0].push_back(t);
        out[1].emplace_back(this->efield(t).real());
        t += this->dt();
        //}
        //else
        //{
            //t = next_pulse_start(t);
        //}
    }

    out[0].shrink_to_fit();
    out[1].shrink_to_fit();
    return out;
}

std::vector<double> PulsetrainParameters::spacing() const
{
    // double tmp;
    // double total = 0;
    std::vector<double> sau( spacing_.size() );
    auto pha = phase();
    for ( size_t i = 0; i < sau.size(); i++ ) {
        sau[i] = spacing_[i] / 0.024188843265; // translate to AU
        // tmp = std::floor( (total + sau[i]) * frequency() / (2 * math::PI) );
        ////std::cout << s_au[i] << ", " << tmp << std::endl;
        // sau[i] = (tmp * 2 * math::PI + pha[i]) / frequency() - total;
        // total += sau[i];
    }
    // however, our spacing isn't exact:
    // divide the spacing by the period (multiply by the frequency?), drop the
    // remainder, and
    // add the phase / 2pi * a period:

    return sau;
}

bool PulsetrainParameters::in_pulse( PetscReal t ) const
{
    //if we are looking at a single pulse, we are always in a pulse
    if ( single_pulse ) return true;
    double tot_time = 0;
    for ( size_t i = 0; i < s_au.size(); i++ ) {
        if ( t < tot_time + s_au[i] && t > pulse_length() + tot_time )
            return false;
        tot_time += s_au[i];
    }
    return true;
}

PetscReal PulsetrainParameters::next_pulse_start( PetscReal t ) const
{
    // PetscReal start = 0;
    PetscReal total = 0;
    for ( size_t i = 0; i < s_au.size(); i++ ) {
        total += s_au[i];
        if ( t < total ) return total;
    }
    return 0;
}

PetscReal PulsetrainParameters::max_time() const
{
    if ( single_pulse ) return pulse_length();
    double tot = 0;
    for ( auto a : s_au )
        tot += a; // one pulse per, the "spacing" is peak to peak.

    return tot + pulse_length();
}

std::vector<double> PulsetrainParameters::phase() const
{
    if ( phase_.size() == spacing_.size() )
        return phase_;
    else
        return std::vector<double>( spacing_.size(), 0 );
}

PetscScalar PulsetrainParameters::efield( PetscReal t ) const
{
    if ( single_pulse ) return LaserParameters::efield( t );
    // find the start of the pulse:
    PetscReal total = 0;
    std::vector<double> ph = phase();
    if ( in_pulse( t ) ) {
        PetscScalar efield = 0;
        for ( size_t i = 0; i < s_au.size(); i++ ) {
            efield += LaserParameters::efield( t - total, ph[i] );
            total += s_au[i];
        }
        efield += LaserParameters::efield( t - total, ph.back() );
        return efield;
    } else
        return 0;
}

void PulsetrainParameters::get_parameters()
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

    if ( opt.isSet( "-pulsetrain_config" ) ) {
        std::string fname;
        opt.get( "-pulsetrain_config" )->getString( fname );
        if ( !opt.importFile( fname.c_str(), '#' ) ) {
            std::cout << "file must exist!" << std::endl;
            throw std::exception();
        }
    }

    if ( opt.isSet( "-pulsetrain_single_pulse" ) )
        single_pulse = true;
    else
        single_pulse = false;

    opt.get( "-pulsetrain_spacing" )->getDoubles( spacing_ );
    opt.get( "-pulsetrain_phase" )->getDoubles( phase_ );
    
    if (opt.isSet( "-pulsetrain_test" ) )
        test_ = true;
    else
        test_ = false;
}

void PulsetrainParameters::save_parameters() const
{
    LaserParameters::save_parameters();
    std::ofstream file;
    file.open( std::string( "./Pulsetrain.config" ) );
    if (single_pulse)
        file << "-pulsetrain_single_pulse" << std::endl;
    file << "-pulsetrain_spacing ";
    for ( size_t i = 0; i < spacing_.size(); i++ ) {
        if ( i == 0 )
            file << spacing_[i];
        else
            file << "," << spacing_[i];
    }
    file << std::endl;
    file << "-pulsetrain_phase ";
    for ( size_t i = 0; i < phase().size(); i++ ) {
        if ( i == 0 )
            file << phase()[i];
        else
            file << "," << phase()[i];
    }
    file << std::endl;
    file.close();
}

std::string PulsetrainParameters::print() const
{
    std::ostringstream out;
    out << LaserParameters::print();
    out << "pulsetrain_spacing: ";
    for ( auto a : spacing_ ) out << a << ",";
    out << std::endl;
    out << "pulsetrain_spacing_au: ";
    for ( auto a : s_au ) out << a << ",";
    out << std::endl;
    out << "pulsetrain_phase: ";
    for ( auto a : phase() ) out << a << ",";
    out << std::endl;
    if (single_pulse)
        out << "pulsetrain_single_pulse" << std::endl;
    return out.str();
}
void PulsetrainParameters::register_parameters()
{
    LaserParameters::register_parameters();

    std::string prefix = "-pulsetrain_";
    opt.overview = "Pulse Train Parameters";
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
    opt.add( "100",
             0,
             1,
             ',',
             "spacing between each pulse (in fs, the last one is repeated)",
             std::string( prefix ).append( "spacing\0" ).c_str() );
    opt.add( "0",
             0,
             1,
             ',',
             "absolute phase of each pulse",
             std::string( prefix ).append( "phase\0" ).c_str() );
    opt.add( "",
             0,
             0,
             0,
             "should we treat this as a 'single pulse'",
             std::string( prefix ).append( "single_pulse\0" ).c_str() );
    opt.add( "./efield.dat",
             0,
             1,
             0,
             "filename to put the laser",
             std::string( prefix ).append( "laser_filename\0" ).c_str() );
    opt.add( "",
             0,
             0,
             0,
             "Test?",
             std::string( prefix ).append( "test\0" ).c_str() );
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
        "Print all inputs and categories for debugging.", // Help description.
        "+d",
        "--debug" // Flag token.
        );
}

