# ebss

The "Energy Basis Schrödinger Solver"

## Setup:

### Homebrew and dependencies

Install [Homebrew](http://brew.sh):

    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"`

Then install needed packages:

    brew install gcc
    brew install openmpi
    brew install openblas
    brew install wget # just to make this a little easier.
    brew install boost --c++11, --with-mpi, --without-single
    brew install gsl

### PETSc and (optional, but easier) SLEPc:

Download and install [PETSc](http://www.mcs.anl.gov/petsc/index.html) in some convenient directory:

    wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.5.3.tar.gz # the lite version doesn't include documentation
    tar zxf petsc-lite-3.5.3.tar.gz
    cd petsc-3.5.3

The following is my configure command for creating the debug version of PETSc. Replace <values in angle brackets> with better values for you:

    ./configure --prefix=<wherever you want to put it> \
                --with-clanguage=cxx \
                --with-c++-support \
                --with-debugging=1 \
                --with-fortran-kernels \
                --with-shared-libraries=true \
                --with-scalar-type=complex \
                --with-blas-lapack-dir=/usr/local/Cellar/openblas/<current_version> \
                --with-mpi=1 \

run the configure, then copy and run each command after it finishes.

#### SLEPc:

Download [SLEPc](http://slepc.upv.es) in some convenient directory:

    wget http://slepc.upv.es/download/download.php?filename=slepc-3.5.3.tar.gz
    mv download.php\?filename=slepc-3.5.3.tar.gz slepc-3.5.3.tar.gz
    tar zxf slepc-3.5.3.tar.gz
    cd slepc-3.5.3

Then the configure command:

    export PETSC_DIR=<prefix from above>; ./configure --prefix=<whereever you want to put it>

run the configure, then copy and run each command after it finishes.

## Modifying and compiling code:

Finally we get to screwing around with ebss itself:

### modify the makefile:

in the root project directory, there is a file: `makefile.include`.  Modify the first two uncommented lines to the prefixes from above, and the third line:

    SLEPC_DIR=<from above for SLEPc>
    PETSC_DIR=<from above for PETSc>
    TDSEBASE=<the full path of the project root directory>

After this, you should be done!
