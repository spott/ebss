#pragma once
#include<boost/function.hpp>
#include<petsc.h>
#include"parameters.hpp"

namespace common
{
   template <typename T>
   Mat populate_matrix(const parameters params,
                       boost::function< bool (int,int) > test,
                       boost::function< T (int,int) > find_value,
                       const unsigned int mat_size,
                       const bool symmetric=true)
   {
      // petsc objects:
      Mat H;

      //Local objects:
      PetscInt start, end;

      MatCreate(params.comm(),&H);
      MatSetType(H, MATMPIAIJ);
      MatSetSizes(H,PETSC_DECIDE,PETSC_DECIDE,
                  mat_size,
                  mat_size);
      MatSetFromOptions(H);

      MatGetOwnershipRange(H, &start, &end);

      T value;

      for (PetscInt i = start; i < end; i++)
         for (PetscInt j = (symmetric ? i : 0u); j < mat_size; j++)

            if (test(i,j))
            {
               value = find_value(i,j);
               MatSetValue(H, i, j, value, INSERT_VALUES);
               if (symmetric)
                  MatSetValue(H, j, i, value, INSERT_VALUES);
            }

      MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY);

      return H;
   }

}
