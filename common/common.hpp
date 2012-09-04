#pragma once

//stl
#include<iostream>
#include<fstream>
#include<vector>
#include<cstdarg>

//boost
#include<boost/function.hpp>

//petsc
#include<petsc.h>

//mine
#include<common/parameters.hpp>

namespace common
{
    template <typename T1, typename T2>
    std::vector<T2> vector_type_change(std::vector<T1> &in)
    {
        std::vector<T2> out(in.size());

        for (size_t i = 0; i < in.size(); i++)
        {
            out[i] = static_cast<T2>(in[i]);
        }

        return out;
    }

    /* Merge Vectors... (and sort?) 
     * sorting requires the > and < operators to be overloaded*/
    template <typename T>
    std::vector<T> merge_vectors(std::vector<T> v1, std::vector<T> v2)
    {
        std::vector<T> out( v1.size() + v2.size() );

        for (size_t i = 0; i < out.size(); i++)
        {
            if (i < v1.size())
                out[i] = v1[i];
            else if (i >= v1.size())
                out[i] = v2[i-v1.size()];
        }
        //C++11 version... would be nice if I could get it to work...
        //auto a = std::move(v1.begin(), v1.end(), out.begin());
        //std::move(v2.begin(), v2.end(), a);

        //delete v1 and v2?
        delete v1;
        delete v2;

        std::sort(out.begin(),out.end());

        return out;
    }

    template <typename T>
    void export_vector_binary(const std::string filename, const std::vector<T> *out)
    {
        std::ios::pos_type size;
        std::ofstream file;
        file.open(filename.c_str(), std::ios::binary | std::ios::out);
        if (file.is_open())
        {
            file.write(reinterpret_cast<const char*>(&(*out)[0]), static_cast<size_t>(sizeof(T)*out->size()));
            file.close();
        }
        else
        {
            std::cerr << "error opening file... does the folder exist?: " << filename << std::endl;
            throw new std::exception();
        }
   };

    template <typename T>
    std::vector<T>* import_vector_binary(const std::string filename)
    {
        //char* buffer;
        std::ios::pos_type size;
        std::ifstream file;
        file.open(filename.c_str(), std::ios::binary | std::ios::in | std::ios::ate);
        std::vector<T> *vec = new std::vector<T>();
        if (file.is_open())
        {
            size = file.tellg();
            if (size != 0)
            {
                std::cerr << size << std::endl;
                vec->resize(size / sizeof(T));
                file.seekg(0, std::ios::beg);

                file.read((char*) vec , size);
                file.close();
            }
            else
            {
                std::cerr << "file is empty!: " << filename << std::endl;
                throw new std::exception();
            }
        }
        else
        {
            std::cerr << "error opening file: " << filename << std::endl;
            throw new std::exception();
        }

        //size_t tsize = sizeof(T);
        //size_t csize = sizeof(char);
        //size_t frac = tsize / csize;

        //for (int i = 0; i < size; i += frac)
            //vec->push_back( static_cast<T>(*(buffer+i)) );

        //delete[] buffer;
        return vec;
   };


    template <typename T>
    std::vector<T>* import_vector_ascii(const std::string filename)
    {
        char* buffer;
        std::vector<T> vec();
        std::ifstream file;
        file.open(filename.c_str(), std::ios::in | std::ios::scientific);
        if (file.is_open())
        {
            T temp;
            while (!file.eof())
            {
                file >> temp;
                vec.push_back(temp);
            }
            file.close();
        }
        else
        {
            std::cerr << "error opening file: " << filename << std::endl;
            throw new std::exception();
        }

        delete[] buffer;
        return vec;
   };

   //template <typename T>
   //Mat populate_matrix(const Parameters params,
                       //boost::function< bool (int,int) > test,
                       //boost::function< T (int,int) > find_value,
                       //const unsigned int mat_size,
                       //const bool symmetric=true)
   //{
      //// petsc objects:
      //Mat H;

      ////Local objects:
      //PetscInt start, end;

      //MatCreate(params.comm(),&H);
      //MatSetType(H, MATMPIAIJ);
      //MatSetSizes(H,PETSC_DECIDE,PETSC_DECIDE,
                  //mat_size,
                  //mat_size);
      //MatSetFromOptions(H);

      //MatGetOwnershipRange(H, &start, &end);

      //T value;

      //for (PetscInt i = start; i < end; i++)
         //for (PetscInt j = (symmetric ? i : 0u); j < mat_size; j++)
            //if (test(i,j))
            //{
               //value = find_value(i,j);
               //MatSetValue(H, i, j, value, INSERT_VALUES);
               //if (symmetric)
                  //MatSetValue(H, j, i, value, INSERT_VALUES);
            //}

      //MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY);
      //MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY);

      //return H;
   //}
}
