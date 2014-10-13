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

class DipoleParameters : public Parameters
{
  public:
    DipoleParameters( int argc, const char** argv, MPI_Comm comm )
        : Parameters( comm )
    {
        register_parameters();

        opt.parse( argc, argv );

        get_parameters();
    }

    DipoleParameters( MPI_Comm comm ) : Parameters( comm )
    {
        // when the parsing will be done elsewhere:
    }

    // call after parsing:
    void get_parameters();

    std::vector<std::string> dipole_filename()
    {
        std::vector<std::string> out;
        out.reserve( 3 * decomp_splits.size() + 1 );
        out.push_back(dipole_filename_);

        std::string df = dipole_filename_.substr( 0, dipole_filename_.find_last_of('.'));

        //This should be enough for now
        const std::vector<char> labels{'a','b','c','d','e','f','g','h','i','j','k'};  

        for (int a = 0; a <= decomp_splits.size(); a++)
          {
            for (int b = a; b <= decomp_splits.size(); b++)
              {
                out.push_back( df + "_" + labels[a] + labels[b] + ".dat");
              }
          }

        return out;
    }
    const std::string after_filename() const
    {
        return after_filename_;
    }
    std::vector<int>& decompositions()
    {
        return decomp_splits;
    }

    PetscScalar find_dipole_moment( Mat& dipole, Vec& psi );

    void find_dipole_moment_decompositions(
        Mat& dipole,
        Vec& psi,
        std::vector<std::ofstream*>& output_vector,
        std::vector<BasisID>& prototype );

    virtual std::string print() const;

    void save_parameters() const;

    PetscReal dt() const
    {
        return dt_ / 0.0241888432652;
    }
    ; // change the fs to A.U.
    PetscReal t_after() const
    {
        return t_after_ / 0.0241888432652;
    }
    ; // ditto

  protected:
    std::vector<int> decomp_splits;
    ez::ezOptionParser opt;
    void register_parameters();
    bool findDipole;
    std::string dipole_filename_;
    std::string after_filename_;
    PetscReal dt_;
    PetscReal t_after_;
    Vec tmp;
};

PetscScalar DipoleParameters::find_dipole_moment( Mat& dipole, Vec& psi )
{
    PetscScalar out;
    if ( tmp == PETSC_NULL ) VecDuplicate( psi, &tmp );
    VecCopy( psi, tmp );
    MatMult( dipole, psi, tmp );
    VecDot( psi, tmp, &out );

    return out;
}

void DipoleParameters::find_dipole_moment_decompositions(
    Mat& dipole,
    Vec& psi,
    std::vector< std::ofstream *>& output_vector,
    std::vector<BasisID>& prototype )
{
    PetscScalar t = this->find_dipole_moment( dipole, psi );
    if (rank() == 0)
	{
        output_vector[0]->write( reinterpret_cast<char*>(&t), sizeof(PetscScalar) );
	}
	MPI_Barrier(MPI_COMM_WORLD);
Vec veca, vecb;
    VecDuplicate( psi, &veca );
    VecDuplicate( psi, &vecb );
    PetscScalar aa, bb, ab;

    std::vector< std::array<int,2> > sections;
    if (decomp_splits.size() != 0)
        for (int i = 0; i <= decomp_splits.size(); i++)
          {
            if (i == decomp_splits.size())
              sections.push_back({decomp_splits.back(), prototype.back().n});
            else if (i == 0)
              sections.push_back({0, decomp_splits[i]});
            else
              sections.push_back({decomp_splits[i-1], decomp_splits[i]});
          }

    int n = 1;
    for ( auto sectiona = sections.begin(); sectiona < sections.end(); sectiona++)
      {
        for ( auto sectionb = sectiona + 1; sectionb < sections.end(); sectionb++)
          {
            // a:
            common::map_function( psi,
                                  [&sectiona, &prototype]( PetscScalar a, PetscInt i ) {
                                    return ( prototype[i].n > (*sectiona)[0] && prototype[i].n <= (*sectiona)[1] ) ? a : 0;
                                  },
                                  veca );
            // b:
            common::map_function( psi,
                                  [&sectionb, &prototype]( PetscScalar a, PetscInt i ) {
                                    return ( prototype[i].n > (*sectionb)[0] && prototype[i].n <= (*sectionb)[1] ) ? a : 0;
                                  },
                                  vecb );

            PetscScalar temp_scalar;

            if ( tmp == PETSC_NULL ) VecDuplicate( psi, &tmp );
            VecCopy( veca, tmp );
            MatMult( dipole, veca, tmp );
            VecDot( vecb, tmp, &temp_scalar );

            if (sectionb == sectiona + 1)
                aa = this->find_dipole_moment( dipole, veca );
			if (sectionb == sections.end() && sectiona == sections.end() - 1)
				bb = this->find_dipole_moment( dipole, vecb );
            ab = temp_scalar;

            //only push back if rank == 0 to preven memory requirement blowup
			if ( rank() == 0) {
				if (sectionb == sectiona + 1)
				{
					output_vector[n]->write( reinterpret_cast<char*>(&aa), sizeof(PetscScalar) );
					n++;
				}
				output_vector[n]->write( reinterpret_cast<char*>(&ab), sizeof(PetscScalar) );
				n++;
				if (sectionb == sections.end() && sectiona == sections.end() - 1)
				{
					output_vector[n]->write( reinterpret_cast<char*>(&bb), sizeof(PetscScalar) );
					n++;
				}
			}
		  }
    }

    VecDestroy( &veca );
    VecDestroy( &vecb );
}

void DipoleParameters::get_parameters()
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

    if ( opt.isSet( "-dipole_config" ) ) {
        std::string fname;
        opt.get( "-dipole_config" )->getString( fname );
        if ( !opt.importFile( fname.c_str(), '#' ) ) {
            std::cout << "file must exist!" << std::endl;
            throw std::exception();
        }
    }
    if ( opt.isSet( "-dipole_filename" ) )
        findDipole = true;
    else
        findDipole = false;

    if ( opt.isSet( "-dipole_decomposition" ) )
        opt.get( "-dipole_decomposition" )->getInts( decomp_splits );

    opt.get( "-dipole_filename" )->getString( dipole_filename_ );
    dipole_filename_ = common::absolute_path( dipole_filename_ );

    opt.get( "-dipole_dt" )->getDouble( dt_ );
    opt.get( "-dipole_t_after" )->getDouble( t_after_ );
    opt.get( "-dipole_after_filename" )->getString( after_filename_ );

    tmp = PETSC_NULL;
}

std::string DipoleParameters::print() const
{
    std::ostringstream out;
    out << "dipole_filename: " << dipole_filename_ << std::endl;
    out << "dipole_after_filename " << after_filename_ << std::endl;
    out << "dipole_dt " << dt_ << std::endl;
    out << "dipole_t_after " << t_after_ << std::endl;
    out << "dipole_decomposition ";
    for( auto i : decomp_splits )
        out << i << ",";
    out << std::endl;
    return out.str();
}
void DipoleParameters::save_parameters() const
{
    std::ofstream file;
    file.open( std::string( "./Dipole.config" ) );
    file << "-dipole_filename " << dipole_filename_ << std::endl;
    file << "-dipole_after_filename " << after_filename_ << std::endl;
    file << "-dipole_dt " << dt_ << std::endl;
    file << "-dipole_t_after " << t_after_ << std::endl;
    file << "-dipole_decomposition ";
    for( auto i : decomp_splits )
        file << i << ",";
    file << std::endl;
    file.close();
}

void DipoleParameters::register_parameters()
{
    std::string prefix = "-dipole_";
    opt.overview = "Dipole Parameters";
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
    opt.add( ".1",
             1,
             1,
             0,
             "dt after (in fs)",
             std::string( prefix ).append( "dt\0" ).c_str() );
    opt.add( "10000",
             1,
             1,
             0,
             "t after (in fs)",
             std::string( prefix ).append( "t_after\0" ).c_str() );
    opt.add( "50",
             1,
             1,
             ',',
             "Split the dipole moment into contributions greater than and less "
             "than x, for each x",
             std::string( prefix ).append( "decomposition\0" ).c_str() );
    opt.add( "Dipole.dat",
             0,
             1,
             0,
             "dipole filename",
             std::string( prefix ).append( "filename\0" ).c_str() );
    opt.add( "after.dat",
             0,
             1,
             0,
             "after the pulses filename for the dipole",
             std::string( prefix ).append( "after_filename\0" ).c_str() );
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
