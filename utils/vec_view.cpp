#include<array>
#include<string>
#include<iostream>
#include<vector>
#include<petsc.h>
#include<common/common.hpp>
#include<common/types.hpp>

int main( int argc, char** argv )
{
    std::array<std::string, 3> fnames;

    if ( argc != 2 && argc != 3)
    {
        std::cout << "one or two arguments only, the vec, or the vec and the prototype" << std::endl;
        std::cout << "only got " << argc << std::endl;
        return 1;
    }
    PetscInitialize( &argc, &argv, PETSC_NULL, PETSC_NULL);

    fnames[0] = std::string( argv[1] );
    if (argc == 3)
        fnames[1] = std::string( argv[2] );

    auto v = common::petsc_binary_read<Vec>(fnames[0], PETSC_COMM_WORLD);
    auto vv = common::Vec_to_vector(v);
    auto proto = common::import_vector_binary<BasisID>(fnames[1]);

    for (auto i = 0; i < vv.size(); i++)
        {
            if (vv[i].real() != 0 or vv[i].imag() != 0)
                std::cout << i << "\t"<< proto[i] << " -- " << vv[i] << "\n";
        }

}

