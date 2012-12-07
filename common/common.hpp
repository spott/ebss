#pragma once

//c stdlib
#include <unistd.h>

//stl
#include<iostream>
#include<fstream>
#include<vector>
#include<cstdarg>
#include<memory>
#include<algorithm>

//petsc
#include<petsc.h>

//mine
#include<common/parameters/Parameters.hpp>

struct BasisID {
    // j = 1 -> 1/2;
    // j = 3 -> 3/2;
    // j != 1,3 -> 0;
    PetscInt n,l,m,j;
    PetscScalar e;
    bool operator<(const BasisID b)
    {
        if (this->l < b.l)
            return true;
        else if (this->l == b.l && this->j < b.j)
            return true;
        else if (this->l == b.l && this->j == b.j && this->n < b.n)
            return true;
        else if (this->l == b.l && this->j == b.j && this->n == b.n && this->m < b.m)
            return true;
        else
            return false;
    }
};
bool operator!=(const BasisID &a, const BasisID &b)
{
    if (a.l != b.l || a.n != b.n || a.m != b.m || a.e != b.e || a.j != b.j)
        return true;
    else
        return false;
}
std::istream& operator>>(std::istream &in, BasisID &b)     //input
{
    PetscReal er, ei;
    in >> b.n >> b.j >> b.l >> b.m >> er >> ei;
    b.e = std::complex<double>(er,ei);
    return in;
}
std::ostream& operator<<(std::ostream &out, const BasisID &b)     //output
{
    out << b.n << ", " << b.j << ", " << b.l << ", " << b.m << ", " << b.e.real() << ", " << b.e.imag();
    return out;
}


namespace common
{

    Vec eigen_balls(Mat mat)
    {
        Vec v;
        MatGetVecs(mat, &v, PETSC_NULL);
        VecSet(v,0.0);
        PetscInt       start = 0, end = 0, row;
        PetscScalar   *array;

        //if (!mat->assembled) SETERRQ(((PetscObject)mat)->comm,PETSC_ERR_ARG_WRONGSTATE,"Not for unassembled matrix");
        //MatCheckPreallocated(mat,1);
        MatGetOwnershipRange(mat, &start, &end);
        VecGetArray(v, &array);
        for(row = start; row < end; ++row) {
            PetscInt           ncols, col;
            const PetscInt    *cols;
            const PetscScalar *vals;

            array[row - start] = 0.0;
            MatGetRow(mat, row, &ncols, &cols, &vals);
            for(col = 0; col < ncols; col++) {
                if (col != row)
                    array[row - start] += std::abs(vals[col]);
            }
            MatRestoreRow(mat, row, &ncols, &cols, &vals);
        }
        VecRestoreArray(v, &array);
        //PetscObjectStateIncrease((PetscObject) v);
        return v;
    }

    std::string absolute_path(std::string rel_path)
    {
        if (rel_path[0] == '.')
        {
            char* a = new char[1025];
            getcwd(a, 1025);
            std::string cwd = std::string(a);
            delete a;
            return cwd.append("/").append(rel_path);
        }
        else
            return rel_path;
    }

    bool file_exists(std::string fname)
    {
        bool ret;
        std::ifstream f(fname);
        if (f.good())
            ret = true;
        else
            ret = false;
        f.close();

        return ret;
    }

    template<typename T>
    T petsc_binary_read(std::string filename, MPI_Comm comm);

    template<>
    Vec petsc_binary_read<Vec>(std::string filename, MPI_Comm comm)
    {
        std::cerr << "importing " << filename << " into vector... " << std::endl;
        if (!file_exists(filename))
        {
            std::cerr << "file doesn't exist" << std::endl;
            throw (std::exception());
        }
        Vec v;
        VecCreate(comm,&v);
        VecSetType(v, VECSTANDARD);
        VecSetFromOptions(v);
        PetscViewer view;
        PetscViewerBinaryOpen(comm, filename.c_str(), FILE_MODE_READ, &view);
        VecLoad(v,view);
        PetscViewerDestroy(&view);
        std::cerr << "done" << std::endl;
        return v;
    }

    template<>
    Mat petsc_binary_read<Mat>(std::string filename, MPI_Comm comm)
    {
        std::cerr << "importing " << filename.c_str() << " into matrix... " << std::endl;
        if (!file_exists(filename))
        {
            std::cerr << "file doesn't exist" << std::endl;
            throw (std::exception());
        }
        Mat v;
        MatCreate(comm,&v);
        MatSetType(v, MATAIJ);
        MatSetFromOptions(v);
        PetscViewer view;
        PetscViewerBinaryOpen(comm, filename.c_str(), FILE_MODE_READ, &view);
        MatLoad(v,view);
        PetscViewerDestroy(&view);
        std::cerr << "done" << std::endl;
        return v;
    }

    void petsc_binary_write(std::string filename, Mat v, MPI_Comm comm)
    {
        PetscViewer view;
        PetscViewerBinaryOpen(comm, filename.c_str(), FILE_MODE_WRITE, &view);
        MatView(v,view);
        PetscViewerDestroy(&view);
    }

    void petsc_binary_write(std::string filename, Vec v, MPI_Comm comm)
    {
        PetscViewer view;
        PetscViewerBinaryOpen(comm, filename.c_str(), FILE_MODE_WRITE, &view);
        VecView(v,view);
        PetscViewerDestroy(&view);
    }

    template <typename T1, typename T2>
    std::vector<T2> vector_type_change(std::vector<T1> *in)
    {
        std::vector<T2> out = new std::vector<T2>(in->size());

        for (size_t i = 0; i < in->size(); i++)
        {
            out[i] = static_cast<T2>((*in)[i]);
        }

        return out;
    }
    template <typename T1, typename T2>
    std::vector<T2> vector_type_change(std::vector<T1> in)
    {
        std::vector<T2> out(in.size());

        for (size_t i = 0; i < in.size(); i++)
        {
            out[i] = static_cast<T2>(in[i]);
        }

        return out;
    }

    template <>
    std::vector<PetscReal> vector_type_change(std::vector<PetscScalar> in)
    {
        std::vector<PetscReal> out(in.size());

        for (size_t i = 0; i < in.size(); i++)
        {
            out[i] = in[i].real();
        }

        return out;
    }
    /* Merge Vectors... (and sort?) 
     * sorting requires the > and < operators to be overloaded*/
    template <typename T>
    std::vector<T> merge_vectors(std::vector<T> v1, std::vector<T> v2)
    {
        std::vector<T> out( v1.size() + v2.size() );

        //C++11 version... would be nice if I could get it to work...
        auto a = std::move(v1.begin(), v1.end(), out.begin());
        std::move(v2.begin(), v2.end(), a);

        //delete v1 and v2?
        delete v1;
        delete v2;

        std::sort(out.begin(),out.end());

        return out;
    }

    template <typename T>
    void export_vector_binary(const std::string filename, const std::vector<T>& out)
    {
        std::ios::pos_type size;
        std::ofstream file;
        file.open(filename.c_str(), std::ios::binary | std::ios::out);
        if (file.is_open())
        {
            file.write(reinterpret_cast<const char*>(&out[0]), static_cast<size_t>(sizeof(T)* out.size()));
            file.close();
        }
        else
        {
            std::cerr << "error opening file... does the folder exist?: " << filename << std::endl;
            throw new std::exception();
        }
    };

    template <typename T>
    std::vector<T> import_vector_binary(const std::string filename)
    {
        std::ios::pos_type size;
        std::ifstream file;
        file.open(filename.c_str(), std::ios::binary | std::ios::in | std::ios::ate);
        std::vector<T> vec;
        if (file.is_open())
        {
            size = file.tellg();
            if (size != 0)
            {
                vec.resize(size / sizeof(T));
                file.seekg(0, std::ios::beg);

                file.read((char*) &vec[0] , size);
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

        return vec;
   };

    template <typename T>
    void export_vector_ascii(const std::string filename, const std::vector<T>& out)
    {
        std::ios::pos_type size;
        std::ofstream file;
        file.open(filename.c_str());
        if (file.is_open())
        {
            file << std::scientific;
            //file.setf(std::ios_base::fixed, std::ios_base::floatfield);
            file.precision(20);
            for (auto a: out)
                file << a << std::endl;

            file.close();
        }
        else
        {
            std::cerr << "file couldn't be opened! " << filename << std::endl;
            throw new std::exception();
        }
    }

    template <typename T>
    std::vector<T> import_vector_ascii(const std::string filename)
    {
        char* buffer;
        std::vector<T> vec;
        std::ifstream file;
        file.open(filename.c_str());
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
   //Vec vector_to_Vec(std::vector< T > vector, MPI_Comm comm)


   //template <>
   //Vec vector_to_Vec<PetscScalar>(std::vector< PetscScalar > vector, MPI_Comm comm)
   //{
       //int size;
       //MPI_Comm_size(comm, &size);

       //Vec v;
       //PetscScalar* array;
       //VecCreate(comm, &v);
       //VecSetSize(v, PETSC_DECIDE, vector.size());
       //VecGetArray(v, array);
       //std::copy(vector.begin(), vector.end(), array);
       //VecRestoreArray(v, array);
       //VecAssemblyBegin(v);
       //VecAssemblyEnd(v);

       //if (size = 1)
       //{
           //Vec v;
           //PetscScalar* array;
           //VecCreate(comm, &v);
           //VecSetSize(v, PETSC_DECIDE, vector.size());
           //VecGetArray(v, array);
           //std::copy(vector.begin(), vector.end(), array);
           //VecRestoreArray(v, array);
           //VecAssemblyBegin(v);
           //VecAssemblyEnd(v);
           //return v;
       //}
       //else
       //{
           //Vec v;
           //PetscScalar* array;
           //VecCreate(comm, &v);
           //VecSetSize(v, PETSC_DECIDE, vector.size());
           //VecGetArray(v, array);

           //VecGet
       //}
   //}


   std::vector<PetscScalar> Vec_to_vector(Vec v)
   {
       MPI_Comm comm;
       int rank;
       int size;
       VecGetSize(v, &size);
       PetscObjectGetComm((PetscObject)v,&comm);
       MPI_Comm_rank(comm, &rank);
       if (rank==0) std::cerr << "entering Vec_to_vector, size = " << size << std::endl;
       Vec seq;
       VecScatter sc;
       VecCreate(comm, &seq);
       //if (rank == 0)
           //VecCreateSeq(comm, size, &seq);
       //if (rank == 0)
            //VecSetSizes(seq, size, size);
       //else
           //VecSetSizes(seq, 0, size);

       VecScatterCreateToZero(v, &sc, &seq);
       if (rank==0) std::cerr << "scattering" << std::endl;
       VecScatterBegin(sc,v,seq,INSERT_VALUES,SCATTER_FORWARD);
       VecScatterEnd(sc,v,seq,INSERT_VALUES,SCATTER_FORWARD);

       if (rank==0) std::cerr << "copying to vector" << std::endl;
       std::vector<PetscScalar> vout(size);
       if (rank==0)
       {
           PetscScalar* a;
           VecGetArray(seq, &a);
           std::copy(a, a+size, vout.begin());
           VecRestoreArray(seq, &a);
       }

       if (rank==0) std::cerr << "destroying vector/scatterer" << std::endl;
       VecScatterDestroy(&sc);
       if (rank == 0) VecDestroy(&seq);
       return vout;
   }


   void printProgBar( double percent )
   {
       percent *= 100;
       std::string bar;

       for(int i = 0; i < 50; i++){
           if( i < (percent/2)){
               bar.replace(i,1,"=");
           }else if( i == (int(percent)/2)){
               bar.replace(i,1,">");
           }else{
               bar.replace(i,1," ");
           }
       }

       std::cout<< "\r" "[" << bar << "] ";
       std::cout.width( 3 );
       std::cout<< percent << "%     " << std::flush;
   };

   template <typename T>
   Mat populate_matrix(const Parameters params,
                       std::function<bool (int,int)> test,
                       std::function<T (int,int)> find_value,
                       const unsigned int mat_size_m,
                       const unsigned int mat_size_n,
                       //const unsigned int diagonal_storage,
                       //const unsigned int offdiag_storage,
                       const bool symmetric=true)
   {
       //if (!symmetric && params.rank() == 0)
           //std::cout << "Calculating for non-symmetric matrix" << std::endl;
       //if (symmetric && params.rank() == 0)
           //std::cout << "Calculating for symmetric matrix" << std::endl;

       // petsc objects:
       Mat H;

       //Local objects:
       PetscInt rowstart, rowend;
       PetscInt colstart, colend;

       MatCreate(params.comm(),&H);
       MatSetType(H, MATMPIAIJ);
       MatSetSizes(H,PETSC_DECIDE,PETSC_DECIDE,
               mat_size_m,
               mat_size_n);

       MatSetUp(H);

       MatGetOwnershipRange(H, &rowstart, &rowend);
       MatGetOwnershipRangeColumn(H, &colstart, &colend);
       PetscInt dnnz[rowend-rowstart];
       PetscInt onnz[rowend-rowstart];
       //find the preallocation functions:
       for (size_t i = rowstart; i < rowend; i++)
       {
           dnnz[i-rowstart] = 0;
           onnz[i-rowstart] = 0;
           for (size_t j = 0; j < mat_size_n; j++)
           {
               if (test(i,j))
               {
                   if (j >= colstart && j < colend)
                       dnnz[i-rowstart]++;
                   else
                       onnz[i-rowstart]++;
               }
           }
       }

       MatMPIAIJSetPreallocation(H, PETSC_NULL, dnnz, PETSC_NULL, onnz);

       MatSetFromOptions(H);

       T value;

       for (PetscInt i = rowstart; i < rowend; i++)
       {
           for (PetscInt j = (symmetric ? i : 0u); j < mat_size_n; j++)
           {
               if (test(i,j))
               {
                   value = find_value(i,j);
                   MatSetValue(H, i, j, value, INSERT_VALUES);
                   if (symmetric)
                       MatSetValue(H, j, i, value, INSERT_VALUES);
               }
           }
           if (params.rank() == 0) 
               printProgBar( double(i - rowstart)/ double(rowend-rowstart) );
       }
       if (params.rank() == 0) std::cout << std::endl;

       MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY);
       MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY);

       return H;
   }

}
