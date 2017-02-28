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

    if ( argc != 2 && argc != 3 ) {
        std::cout << "one or two arguments only, the matrix, or the matrix and "
                     "the prototype"
                  << std::endl;
        std::cout << "only got " << argc << std::endl;
        return 1;
    }
    PetscInitialize( &argc, &argv, PETSC_NULL, PETSC_NULL );

    fnames[0]                  = std::string( argv[1] );
    if ( argc == 3 ) fnames[1] = std::string( argv[2] );

    Mat mat1;
    MatCreate( PETSC_COMM_WORLD, &mat1 );
    MatSetType( mat1, MATAIJ );
    MatSetFromOptions( mat1 );


    PetscViewer view;
    PetscViewerBinaryOpen( PETSC_COMM_WORLD, fnames[0].c_str(), FILE_MODE_READ,
                           &view );
    MatLoad( mat1, view );

    Vec temp;
    Vec result;
    MatGetVecs( mat1, &temp, &result );

    // std::vector<BasisID> proto;
    // if (argc == 3)
    auto proto = common::import_vector_binary<BasisID>( fnames[1] );
    // else
    //     {
    //         PetscInt s;
    //         VecGetSize(result, &s);
    //         for (auto i = 0; i < s; i++)
    //             {
    //                 proto.push_back(BasisID(i,0,0,0,std::complex<double>(0,0)));
    //             }
    //     }
    while ( true ) {
        VecZeroEntries( temp );
        int column = 0;
        std::cout << "which column? ";
        std::cin >> column;
        VecSetValue( temp, column, 1, INSERT_VALUES );
        VecAssemblyBegin( temp );
        VecAssemblyEnd( temp );
        MatMult( mat1, temp, result );
        auto v = common::Vec_to_vector( result );
        for ( int i = 0; i < v.size(); i++ ) {
            if ( v[i].real() != 0.0 or v[i].imag() != 0.0 )
                std::cout << i << "\t" << proto[i] << " -- " << v[i] << "\n";
        }
        // VecView(result, PETSC_VIEWER_STDOUT_WORLD);
    }
}
