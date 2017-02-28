#include <array>
#include <common/common.hpp>
#include <common/types.hpp>
#include <iostream>
#include <petsc.h>
#include <string>
#include <vector>

int main( int argc, char** argv )
{
    std::array<std::string, 3> fnames;

    if ( argc != 4 && argc != 3 ) {
        std::cout << "two or three arguments expected: filename for matrix 1, "
                     "matrix 2, and matrix 2 - matrix 1 "
                  << std::endl;
        std::cout << "only got " << argc << std::endl;
        return 1;
    }
    PetscInitialize( &argc, &argv, PETSC_NULL, PETSC_NULL );

    fnames[0] = std::string( argv[1] );
    fnames[1] = std::string( argv[2] );
    fnames[2] = std::string( argv[3] );

    Mat mat1;
    MatCreate( PETSC_COMM_WORLD, &mat1 );
    MatSetType( mat1, MATAIJ );
    MatSetFromOptions( mat1 );
    Mat mat2;
    MatCreate( PETSC_COMM_WORLD, &mat2 );
    MatSetType( mat2, MATAIJ );
    MatSetFromOptions( mat2 );


    PetscViewer view;
    PetscViewerBinaryOpen( PETSC_COMM_WORLD, fnames[0].c_str(), FILE_MODE_READ,
                           &view );
    MatLoad( mat1, view );
    PetscViewerBinaryOpen( PETSC_COMM_WORLD, fnames[1].c_str(), FILE_MODE_READ,
                           &view );
    MatLoad( mat2, view );

    MatAYPX( mat1, -1, mat2, SAME_NONZERO_PATTERN );

    int row, cols;
    MatGetSize( mat1, &row, &cols );
    Vec v;
    MatGetVecs( mat1, &v, PETSC_NULL );
    std::vector<PetscInt> idx( cols );
    // PetscInt* idx = new PetscInt[col];
    MatGetRowMaxAbs( mat1, v, idx.data() );

    // VecView(v, PETSC_VIEWER_STDOUT_WORLD);
    auto proto = common::import_vector_binary<BasisID>(
        "500/120000/300_hamiltonian/vector_prototype.dat" );

    PetscScalar* p;
    VecGetArray( v, &p );

    for ( size_t i = 0; i < proto.size(); ++i ) {
        std::cout << proto[i] << " <--> " << proto[idx[i]] << " = " << p[i]
                  << std::endl;
    }
    // for(auto a: idx)
    //{
    // std::cout << a << std::endl;
    //}

    PetscViewerBinaryOpen( PETSC_COMM_WORLD, fnames[2].c_str(), FILE_MODE_WRITE,
                           &view );
    MatView( mat1, view );
}
