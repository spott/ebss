
#include<petsc.h>

//stl:
#include<complex>
#include<vector>
#include<cassert>

//mine:
#include<common/parameters/HamiltonianParameters.hpp>
#include<common/parameters/NonlinearParameters.hpp>
#include<common/parameters/StateParameters.hpp>
#include<common/math.hpp>
#include<common/common.hpp>


Vec psi( int order, std::vector< double >::const_iterator frequencies_begin, std::vector<double>::const_iterator frequencies_end, PetscScalar wg, Vec& H0, Mat& D, Vec& psi0, Vec& mask);
Vec psi_conjugate( int order, std::vector< double >::const_iterator frequencies_begin, std::vector<double>::const_iterator frequencies_end, PetscScalar wg, Vec& H0, Mat& D, Vec& psi0, Vec& mask);

int main( int argc, const char** argv )
{
    int ac = argc;
    char** av = new char*[argc];
    for (size_t i = 0; i < argc; i++)
    {
        av[i] = new char[strlen(argv[i])+1];
        std::copy(argv[i], argv[i] + strlen(argv[i])+1, av[i]);
    }
    PetscInitialize(&ac, &av, PETSC_NULL, PETSC_NULL);

    //the parameters (where to find the hamiltonian)
    NonlinearParameters nparams(argc, argv, MPI_COMM_WORLD);
    HamiltonianParameters<double> params(argc, argv, MPI_COMM_WORLD);
    StateParameters sparams(argc, argv, MPI_COMM_WORLD);
    //the state: for removing states (such as when the ground state isn't the ground state calculated.

    //the parameters for this:

    std::cout << nparams.print();


    //read in the matrices
    Mat D = params.read_dipole_matrix();
    Vec H0 = params.read_energy_eigenvalues();

    //create a mask for bound states:
    Vec mask = common::map_function(H0, [](PetscScalar in) { return 1; });
    //int start, end;
    //PetscScalar* a;
    //VecGetOwnershipRange(mask, &start, &end);
    //VecGetArray(mask, &a);
    //for(int i = start; i < end; ++i)
        //if (a[i - start] == 0.0)
            //MatZeroRowsColumns(
    //View the hamiltonian
    //VecView(H0, PETSC_VIEWER_STDOUT_WORLD);

    VecShift(H0, std::complex<double>(0,-.00001));
    Vec psi0;
    VecDuplicate(H0, &psi0);
    VecSetValue(psi0, 0, 1., INSERT_VALUES);
    VecAssemblyBegin(psi0);
    VecAssemblyEnd(psi0);

    //find wg (the ground state energy)
    PetscScalar wg;
    VecDot(psi0, H0, &wg);
    if (params.rank() == 0) std::cout << "wg: " << wg << std::endl;

    //the size of the frequency run:
    const int num_freqs = 200;
    const double freq_step = 0.005;

    // get the frequencies list from nonlinear params:
    auto freqs = nparams.freqs();

    //create the vectors that will store the chi's
    //std::vector< std::complex<double> > chi1;
    //chi1.reserve(frequencies.size());
    //std::vector< std::complex<double> > chi3;
    //chi3.reserve(frequencies.size());
    //std::vector< std::complex<double> > chi5;
    //chi5.reserve(frequencies.size());

    std::vector< std::vector< std::complex<double> > > chi1_data( nparams.chi1s().size() );
    for( auto a : chi1_data) a.reserve( freqs.size() );
    std::vector< std::vector< std::complex<double> > > chi3_data( nparams.chi3s().size() );
    for( auto a : chi3_data) a.reserve( freqs.size() );
    std::vector< std::vector< std::complex<double> > > chi5_data( nparams.chi5s().size() );
    for( auto a : chi5_data) a.reserve( freqs.size() );

    for( int f = 0; f < freqs.size(); ++f)
    {
        PetscScalar t1,t2,t3,t4,t5,t6;
        Vec c;
        VecDuplicate(H0, &c);
        //we need to loop through all permutations of all chi types:

        //first we need to generate the correct vectors:
        //Chi1
        for (auto i = nparams.chi1s().begin(); i != nparams.chi1s().end(); ++i)
        {
            std::vector< double > freq{(*i) * freqs[f]};
            Vec p0  = psi(0, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
            Vec p1  = psi(1, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
            Vec p1c = psi_conjugate(1, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
            Vec p0c = psi_conjugate(0, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);

            // chi1 = <\psi^(0) | D | \psi^(1)>
            MatMult(D, p1, c);
            VecDot(p0c, c, &t1);
            // chi1 = <\psi^(1) | D | \psi^(0)>
            MatMult(D, p0, c);
            VecDot(p1c, c , &t2);
            chi1_data[i - nparams.chi1s().begin()].push_back((t1+t2));
            VecDestroy(&p1);
            VecDestroy(&p0);
            VecDestroy(&p1c);
            VecDestroy(&p0c);
        }

        //Chi3
        for (auto i = nparams.chi3s().begin(); i != nparams.chi3s().end(); ++i)
        {
            std::sort((*i).begin(), (*i).end());
            size_t multiplicity = 1;
            std::array<int, 3> ts{0,0,0};
            for(auto m : (*i))
            {
                if (m == -1)
                    ts[0]++;
                if (m == 1)
                    ts[2]++;
                if (m == 0)
                    ts[1]++;
            }
            for (auto m: ts)
                multiplicity *= math::factorial(m);

            PetscScalar result;
            do {

                std::vector< double > freq{(*i)[0] * freqs[f], (*i)[1] * freqs[f], (*i)[2] * freqs[f]};
                Vec p3  = psi(3, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
                Vec p2  = psi(2, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
                Vec p0  = psi(0, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
                Vec p1  = psi(1, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
                Vec p3c = psi_conjugate(3, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
                Vec p2c = psi_conjugate(2, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
                Vec p1c = psi_conjugate(1, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
                Vec p0c = psi_conjugate(0, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);

                // chi3 = <\psi^(0) | D | \psi^(3)>
                MatMult(D, p3, c);
                VecDot(p0c ,c,  &t1);
                //      + <\psi^(1) | D | \psi^(2)>
                MatMult(D, p2, c);
                VecDot(p1c ,c,  &t2);
                //      + <\psi^(2) | D | \psi^(1)>
                MatMult(D, p1, c);
                VecDot(p2c ,c,  &t3);
                ////      + <\psi^(3) | D | \psi^(0)>
                MatMult(D, p0, c);
                VecDot(p3c ,c,  &t4);
                result += (t1 + t2 + t3 + t4);
                VecDestroy(&p3);
                VecDestroy(&p2);
                VecDestroy(&p1);
                VecDestroy(&p0);
                VecDestroy(&p3c);
                VecDestroy(&p2c);
                VecDestroy(&p1c);
                VecDestroy(&p0c);
            } while (std::next_permutation( (*i).begin(), (*i).end() ) );
            chi3_data[i-nparams.chi3s().begin()].push_back(result * static_cast<double>(multiplicity) / static_cast<double>(math::factorial(3)));
        }
        
        //Chi5
        for (auto i = nparams.chi5s().begin(); i != nparams.chi5s().end(); ++i)
        {
            std::sort((*i).begin(), (*i).end());
            size_t multiplicity = 1;
            std::array<int, 3> ts{0,0,0};
            for(auto m : (*i))
            {
                if (m == -1)
                    ts[0]++;
                if (m == 0)
                    ts[1]++;
                if (m == 1)
                    ts[2]++;
            }
            for (auto m: ts)
                multiplicity *= math::factorial(m);
            PetscScalar result;
            do {
                std::vector< double > freq{(*i)[0] * freqs[f], (*i)[1] * freqs[f], (*i)[2] * freqs[f], (*i)[3] * freqs[f], (*i)[4] * freqs[f]};
                Vec p5  = psi(5, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
                Vec p4  = psi(4, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
                Vec p3  = psi(3, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
                Vec p2  = psi(2, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
                Vec p0  = psi(0, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
                Vec p1  = psi(1, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
                Vec p5c = psi_conjugate(5, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
                Vec p4c = psi_conjugate(4, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
                Vec p3c = psi_conjugate(3, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
                Vec p2c = psi_conjugate(2, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
                Vec p1c = psi_conjugate(1, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
                Vec p0c = psi_conjugate(0, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);

                // chi5 = <\psi^(0) | D | \psi^(5)>
                MatMult(D, p5, c);
                VecDot(p0c ,c,  &t1);
                //      + <\psi^(1) | D | \psi^(4)>
                MatMult(D, p4, c);
                VecDot(p1c ,c,  &t2);
                //      + <\psi^(2) | D | \psi^(3)>
                MatMult(D, p3, c);
                VecDot(p2c ,c,  &t3);
                //      + <\psi^(2) | D | \psi^(3)>
                MatMult(D, p2, c);
                VecDot(p3c ,c,  &t4);
                //      + <\psi^(4) | D | \psi^(1)>
                MatMult(D, p1, c);
                VecDot(p4c ,c,  &t5);
                //      + <\psi^(5) | D | \psi^(0)>
                MatMult(D, p0, c);
                VecDot(p5c ,c,  &t6);
                result += (t1 + t2 + t3 + t4);
                VecDestroy(&p5);
                VecDestroy(&p4);
                VecDestroy(&p3);
                VecDestroy(&p2);
                VecDestroy(&p1);
                VecDestroy(&p0);
                VecDestroy(&p5c);
                VecDestroy(&p4c);
                VecDestroy(&p3c);
                VecDestroy(&p2c);
                VecDestroy(&p1c);
                VecDestroy(&p0c);
            } while (std::next_permutation( (*i).begin(), (*i).end() ) );
            chi5_data[i-nparams.chi5s().begin()].push_back( result * static_cast<double>(multiplicity) / static_cast<double>(math::factorial(5)));
        }

        if (params.rank() == 0) std::cout << "..." << f << "(" << freqs[f] << ")" << std::flush;
        VecDestroy(&c);
    }

    //export vectors.
    std::stringstream ss;
    for( auto a = chi1_data.begin(); a != chi1_data.end(); ++a)
    {
        ss << "chi1_" << nparams.chi1s()[a - chi1_data.begin()] << ".dat";
        common::export_vector_binary( ss.str(), *a);
        ss.str() = "";
    }
    ss.str("");

    for( auto a = chi3_data.begin(); a != chi3_data.end(); ++a)
    {
        ss << "chi3_" ;
        for( auto b : nparams.chi3s()[a - chi3_data.begin()])
            ss<< b << "_";
        ss << ".dat";
        common::export_vector_binary( ss.str(), *a);
        ss.str() = "";
    }

    ss.str("");
    for( auto a = chi5_data.begin(); a != chi5_data.end(); ++a)
    {
        ss << "chi5_" ;
        for( auto b : nparams.chi5s()[a - chi5_data.begin()])
            ss << b << "_";
        ss << ".dat";
        common::export_vector_binary( ss.str(), *a);
        ss.str() = "";
    }
    //for( auto a = chi3_data.begin(); a != chi3_data.end(); ++a)
        //common::export_vector_binary( std::string("chi1_") + std::string(nparams.chi1s()[a - chi1_data.begin()]) + std::string(".dat") , a);
    //common::export_vector_binary( "chi3.dat" , chi3);
    //common::export_vector_binary( "chi5.dat" , chi5);

    std::cout << std::endl << "done with general code." <<std::endl;

    /*****************
     * Chi1!
     *****************/

    //temporary denominator vector
    //Vec d;
    //VecDuplicate(H0, &d);

    ////tmp vectors for matrix multiplication:
    //Vec tmp, tmp2;
    //VecDuplicate(H0, &tmp);
    //VecDuplicate(H0, &tmp2);

    //// the output vector for chi1 "by the book"
    //std::vector< std::complex<double> > out;
    //out.reserve(10000);

    //for (int i = 0; i < 10000; ++i)
    //{
        ////calculate chi1 "by the book"
        //double freq = .0001*i;

        //PetscScalar t1, t2;

        //VecCopy(H0, d);
        //VecShift(d, -wg);
        //VecShift(d, -freq);
        //VecReciprocal(d);

        //MatMult(D, psi0, tmp);
        //VecPointwiseMult(tmp2, d, tmp);
        //MatMult(D, tmp2, tmp);
        //VecDot(psi0, tmp, &t1);


        //VecCopy(H0, d);
        //VecShift(d, -wg);
        //VecShift(d, freq);
        //VecConjugate(d);
        //VecReciprocal(d);

        //MatMult(D, psi0, tmp);
        //VecPointwiseMult(tmp2, d, tmp);
        //MatMult(D, tmp2, tmp);
        //VecDot(psi0, tmp, &t2);

        //out.push_back(pre*(t1+t2));
    //}
    //common::export_vector_binary( "chi1_book.dat" ,out);

    /*****************
     * Chi3!
     *****************/

    ////three denominators:
    //Vec d1,d2,d3;
    //VecDuplicate(H0, &d1);
    //VecDuplicate(H0, &d2);
    //VecDuplicate(H0, &d3);

    ////we will reuse the vector for chi1
    //out.clear();

    //for( int i = 0; i < 10000; ++i)
    //{
        //double freq = .0001*i;

        ////four terms, each with 3 denominators
        //PetscScalar t1,t2,t3,t4;

        ////term 1:
        //VecCopy(H0, d1);
        //VecShift(d1, -wg);
        //VecShift(d1, -freq);

        //VecCopy(d1, d2);
        //VecShift(d2, -freq);

        //VecCopy(d2, d3);
        //VecShift(d3, -freq);

        //VecReciprocal(d1);
        //VecReciprocal(d2);
        //VecReciprocal(d3);

        //MatMult(D, psi0, tmp);
        //VecPointwiseMult(tmp, d1, tmp);
        //MatMult(D, tmp, tmp2);
        //VecPointwiseMult(tmp, d2, tmp2);
        //MatMult(D, tmp, tmp2);
        //VecPointwiseMult(tmp, d3, tmp2);
        //MatMult(D, tmp, tmp2);

        //VecDot(psi0, tmp2 , &t1);

        ////term 2:
        //VecCopy(H0, d1);
        //VecShift(d1, -wg);
        //VecShift(d1, -freq);

        //VecCopy(d1, d2);
        //VecShift(d2, -freq);

        //VecCopy(H0, d3);
        //VecShift(d3, -wg);
        //VecShift(d3, freq);
        //VecConjugate(d3);

        //VecReciprocal(d1);
        //VecReciprocal(d2);
        //VecReciprocal(d3);

        //MatMult(D, psi0, tmp);
        //VecPointwiseMult(tmp, d1, tmp);
        //MatMult(D, tmp, tmp2);
        //VecPointwiseMult(tmp, d2, tmp2);
        //MatMult(D, tmp, tmp2);
        //VecPointwiseMult(tmp, d3, tmp2);
        //MatMult(D, tmp, tmp2);

        //VecDot(psi0, tmp2 , &t2);

        ////term 3:
        //VecCopy(H0, d1);
        //VecShift(d1, -wg);
        //VecShift(d1, -freq);

        //VecCopy(H0, d2);
        //VecShift(d2, -wg);
        //VecShift(d2, freq);
        //VecShift(d2, freq);
        //VecConjugate(d2);

        //VecCopy(H0, d3);
        //VecShift(d3, -wg);
        //VecShift(d3, freq);
        //VecConjugate(d3);

        //VecReciprocal(d1);
        //VecReciprocal(d2);
        //VecReciprocal(d3);

        //MatMult(D, psi0, tmp);
        //VecPointwiseMult(tmp, d1, tmp);
        //MatMult(D, tmp, tmp2);
        //VecPointwiseMult(tmp, d2, tmp2);
        //MatMult(D, tmp, tmp2);
        //VecPointwiseMult(tmp, d3, tmp2);
        //MatMult(D, tmp, tmp2);

        //VecDot(psi0, tmp2 , &t3);

        ////term 4:
        //VecCopy(H0, d1);
        //VecShift(d1, -wg);
        //VecShift(d1, freq);
        //VecShift(d1, freq);
        //VecShift(d1, freq);
        //VecConjugate(d1);

        //VecCopy(H0, d2);
        //VecShift(d2, -wg);
        //VecShift(d2, freq);
        //VecShift(d2, freq);
        //VecConjugate(d2);

        //VecCopy(H0, d3);
        //VecShift(d3, -wg);
        //VecShift(d3, freq);
        //VecConjugate(d3);

        //VecReciprocal(d1);
        //VecReciprocal(d2);
        //VecReciprocal(d3);

        //MatMult(D, psi0, tmp);
        //VecPointwiseMult(tmp, d1, tmp);
        //MatMult(D, tmp, tmp2);
        //VecPointwiseMult(tmp, d2, tmp2);
        //MatMult(D, tmp, tmp2);
        //VecPointwiseMult(tmp, d3, tmp2);
        //MatMult(D, tmp, tmp2);

        //VecDot(psi0, tmp2 , &t4);

        //if (params.rank() == 0) std::cout << "..." << i << std::flush;
        //out.push_back(pre*(t1+t2+t3+t4));
    //}

    //common::export_vector_binary( "chi3_book.dat" ,out);


    PetscFinalize();
}

// Perturbative calculations:  \psi^(n) / E(w) etc.
// Each version takes "n" frequencies in a vector of real values

//template< typename Iterator >
Vec psi( int order, std::vector<double>::const_iterator frequencies_begin, std::vector<double>::const_iterator frequencies_end, PetscScalar wg, Vec& H0, Mat& D, Vec& psi0, Vec& mask)
{
    //we must have enough frequencies to do the calculation:
    assert( frequencies_end - frequencies_begin >= order );

    //the out vector (the one we return)
    Vec out;

    //a temporary vector
    Vec tmp;

    //make sure that out and tmp have the correct memory layout.
    VecDuplicate(psi0, &out);
    VecDuplicate(H0, &tmp);

    //out starts as psi0
    VecCopy(psi0, out);

    //do this "order" times
    for (int i = 1; i <= order; ++i)
    {
        // D \psi0 = tmp
        MatMult(D, out, tmp);
        // \psi0 = tmp
        VecCopy(tmp, out);
        // tmp = H0
        VecCopy(H0, tmp);
        // tmp = tmp - wg
        VecShift(tmp, -wg);

        //subtract the laser frequencies from H0 - wg
        for (auto a = frequencies_begin + order - i; a < frequencies_begin + order; ++a)
        {
            // tmp = tmp - \sum_a \omega_a
            VecShift(tmp, -(*a));
        }
        // tmp = 1/tmp
        VecReciprocal(tmp);
        // out = tmp * out
        VecPointwiseMult(out, tmp, out);
        VecPointwiseMult(out, mask, out);
    }
    
    //destroy the temporary
    VecDestroy(&tmp);

    return out;
}

Vec psi_conjugate( int order, std::vector< double >::const_iterator frequencies_begin, std::vector<double>::const_iterator frequencies_end, PetscScalar wg, Vec& H0, Mat& D, Vec& psi0, Vec& mask)
{
    assert( frequencies_end - frequencies_begin >= order );
    Vec out;
    Vec tmp;
    VecDuplicate(psi0, &out);
    VecDuplicate(H0, &tmp);
    VecCopy(psi0, out);

    for (int i = 1; i <= order; ++i)
    {
        MatMult(D, out, tmp);
        VecCopy(tmp, out);
        VecCopy(H0, tmp);
        VecShift(tmp, -wg);
        for (auto a = frequencies_begin + i - 1; a >= frequencies_begin; --a)
        {
            VecShift(tmp, (*a) );
        }
        VecReciprocal(tmp);
        VecPointwiseMult(out, tmp, out);
        VecPointwiseMult(out, mask, out);
    }

    VecDestroy(&tmp);

    return out;
}
