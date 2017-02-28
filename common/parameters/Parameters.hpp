#pragma once


#include <common/ezOptionParser.hpp>

#include <algorithm>
#include <iostream>
#include <petsc.h>
#include <string>

class Parameters
{
  public:
    Parameters( MPI_Comm comm ) : comm_( comm ){};

    MPI_Comm comm() const { return comm_; };

    int rank() const;
    int size() const;

  protected:
    // virtual void init_from_file();
    // virtual void save_to_file();
    MPI_Comm comm_;
    // std::string params_filename;
};

int Parameters::rank() const
{
    int i = -1;
    MPI_Comm_rank( this->comm_, &i );

    if ( i < 0 ) throw( std::exception() );
    return i;
}

int Parameters::size() const
{
    int i = -1;
    MPI_Comm_size( this->comm_, &i );

    if ( i < 0 ) throw( std::exception() );
    return i;
}
