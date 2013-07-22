
#include<petsc.h>

//stl:
#include<complex>
#include<vector>
#include<cassert>

//mine:
#include<common/parameters/HamiltonianParameters.hpp>
#include<common/math.hpp>
#include<common/common.hpp>


//template< typename Iterator >
//Vec psit( int order, Iterator frequencies_begin, Iterator frequencies_end, PetscScalar wg, Vec& H0, Mat& D, Vec& psi0)
//{
    //assert( std::abs(frequencies_end - frequencies_begin) >= order );
    //Vec out;
    //Vec tmp;
    //VecDuplicate(psi0, &out);
    //VecDuplicate(H0, &tmp);
    //VecCopy(psi0, out);

    //for (int i = 1; i <= order; ++i)
    //{
        //MatMult(D, out, tmp);
        //VecCopy(tmp, out);
        //VecCopy(H0, tmp);
        //VecShift(tmp, wg);
        //for (auto a = frequencies_begin + order - i; a < frequencies_begin + order + 1; ++a)
        //{
            //VecShift(tmp, *a);
        //}
        //VecReciprocal(tmp);
        //VecPointwiseMult(out, tmp, out);
    //}
    
    //VecDestroy(&tmp);

    //return out;
//}
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
    HamiltonianParameters<double> params(argc, argv, MPI_COMM_WORLD);

    //the state: for removing states (such as when the ground state isn't the ground state calculated.
    //StateParameters sparams(argc, argv, MPIC_COMM_WORLD);

    //the parameters for this:



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

    //pre factor for "Chi" vs polarizability numbers.  
    double N = 3.98172e-6;
    double pre = N * math::PI;
    //currently set to 1, so I'm calculating polarizabilities
    pre = 1.;

    //find wg (the ground state energy)
    PetscScalar wg;
    VecDot(psi0, H0, &wg);

    if (params.rank() == 0) std::cout << "wg: " << wg << std::endl;

    //the size of the frequency run:
    const int num_freqs = 200;
    const double freq_step = 0.005;

    // create the frequencies list:
    std::vector< units::Meter > nm{ 800._nm, 400._nm, 263._nm, 266.666_nm, 1200._nm, 1800._nm, 4000._nm};
    std::vector< double > frequencies{ 0. };

    for( auto a : nm )
        frequencies.push_back( units::toEnergy( a ) );

    for( auto a : frequencies )
        std::cout << a << ", ";
    std::cout << std::endl;

    //create the vectors that will store the chi's
    std::vector< std::complex<double> > chi1;
    chi1.reserve(frequencies.size());
    std::vector< std::complex<double> > chi3;
    chi3.reserve(frequencies.size());
    std::vector< std::complex<double> > chi5;
    chi5.reserve(frequencies.size());
    for( int i = 0; i < frequencies.size(); ++i)
    {
        //the frequency list, currently, all frequencies are the same (for 3rd harmonic generation), 
        //later I will need this to look at the Kerr effect and others
        std::vector<double> freq{frequencies[i], frequencies[i], frequencies[i], frequencies[i], frequencies[i]};
        PetscScalar t1,t2,t3,t4,t5,t6;
        Vec c;
        VecDuplicate(H0, &c);

        //perterbative vectors.  the only difference with psi_conjugate is that it takes the positive of the frequencies 
        //(this is necessary so that every term has the same exponential phase term, so that said term can be ignored). 
        //(the actual conjugation happens in VecDot).
        Vec p5  = psi(5, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
        Vec p4  = psi(4, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
        Vec p3  = psi(3, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
        Vec p2  = psi(2, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
        Vec p1  = psi(1, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
        Vec p0  = psi(0, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
        Vec p5c = psi_conjugate(5, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
        Vec p4c = psi_conjugate(4, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
        Vec p3c = psi_conjugate(3, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
        Vec p2c = psi_conjugate(2, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
        Vec p1c = psi_conjugate(1, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);
        Vec p0c = psi_conjugate(0, freq.cbegin(), freq.cend(), wg, H0, D, psi0, mask);

        // chi1 = <\psi^(0) | D | \psi^(1)>
        MatMult(D, p1, c);
        VecDot(p0c, c, &t1);
        // chi1 = <\psi^(1) | D | \psi^(0)>
        MatMult(D, p0, c);
        VecDot(p1c, c , &t2);
        chi1.push_back(pre*(t1+t2));

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
        chi3.push_back(pre*(t1+t2+t3+t4));

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
        chi5.push_back(pre*(t1+t2+t3+t4+t5+t6));


        if (params.rank() == 0) std::cout << "..." << i << std::flush;
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
        VecDestroy(&c);
    }

    //export vectors.
    common::export_vector_binary( "chi1.dat" , chi1);
    common::export_vector_binary( "chi3.dat" , chi3);
    common::export_vector_binary( "chi5.dat" , chi5);

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
