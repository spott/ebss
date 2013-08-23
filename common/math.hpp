#pragma once
//ebss
#include<common/parameters/Parameters.hpp>
#include<common/common.hpp>
#include<common/special/coulomb/complex_functions.H>
#include<common/special/coulomb/cwfcomp.cpp>
#include<common/types.hpp>
//#include<common/special/coulomb/test_rec_rel.>

//fftw3:
extern "C" {
#include <fftw3.h>
}

//stl:
#include<complex>
#include<vector>
#include<cmath>
#include<cassert>
#include<tuple>

//gsl:
#include<gsl/gsl_sf_coupling.h>
#include<gsl/gsl_sf_coulomb.h>
#include<gsl/gsl_sf_legendre.h>



namespace math{
//Constants:

constexpr PetscReal PI = std::atan(1.0)*4.0;
const PetscReal C = 137.035999;

template <typename T> inline constexpr
int signum(T x, std::false_type is_signed) {
    return T(0) < x;
}

template <typename T> inline constexpr
int signum(T x, std::true_type is_signed) {
    return (T(0) < x) - (x < T(0));
}

template <typename T> inline constexpr
int signum(T x) {
    return signum(x, std::is_signed<T>());
}

template <typename scalar>
scalar CGCoefficient(const BasisID &init, const BasisID &fin)
{
    scalar out = gsl_sf_coupling_3j(init.l*2, 2, fin.l*2, 0, 0, 0);
    out *= out;
	out *= std::sqrt(( 2 * init.l +1 ) * ( 2 * fin.l +1 ) );
    return out;
}

//interpolation class: input an iterator

//template <size_t order>
class interpolate
{
    public:
    interpolate( const std::vector< double > & y, const std::vector< double > & x ) : y_(y), x_(x), coefficient(x.size()) {
        //we need to set up the second derivative vector?
        coefficient.front() = 3 * ( (y_[1] - y_[0])/(std::pow(x_[1] - x_.front(),2) ) ); //natural boundary conditions.
        coefficient.back() = 3 * ( (y_.back() - y_[y_.size() - 2])/(std::pow(x_.back() - x_[x_.size() - 2],2) ) ); //natural boundary conditions.
        for( size_t i = 1; i < coefficient.size()-1; ++i)
            coefficient[i] = 3 * ( (y_[i] - y_[i-1])/(std::pow(x_[i] - x_[i-1],2)) + (y_[i+1] - y_[i])/(std::pow(x_[i+1] - x_[i],2)) );

        //populate and solve the tridiagonal system, saving the solution in coefficient_?
        std::vector<double> u(x_.size()); //upper
        std::vector<double> d(x_.size()); //diagonal
        std::vector<double> l(x_.size()); //lower == upper
        d[0] = 2. / (x_[1] - x_[0]);
        u[0] = 1. / (x_[1] - x_[0]);
        l[0] = 1. / (x_[1] - x_[0]);

        for( size_t i = 1; i < x_.size() - 1; ++i)
        {
            u[i] = 1. / (x_[i+1] - x_[i]);
            l[i] = 1. / (x_[i+1] - x_[i]);
            d[i] = 2. * (u[i] + u[i-1]);
        }
        d.back() = 2. * u[x_.size() - 2];

        for( size_t i = 1; i < x_.size(); ++i)
        {
            double m = 1.0 / (d[i] - l[i-1] * u[i - 1]);
            u[i] = u[i] * m;
            coefficient[i] = (coefficient[i] - l[i-1] * coefficient[i - 1]) * m;
        }

        for (size_t i = x_.size() - 1; i --> 0; )
            coefficient[i] = coefficient[i] - u[i] * coefficient[i + 1];

        //now, coefficient is a list of coefficients that I can use.
    };


    double operator() ( double x ) {
        size_t p = find_nearest_point(x);

        if (p == x_.size())
            return y_.back();

        auto t = (x - x_[p]) / (x_[p+1] - x_[p]);
        auto a = coefficient[p] * (x_[p+1] - x_[p]) - (y_[p+1] - y_[p]);
        auto b = -coefficient[p+1] * (x_[p+1] - x_[p]) + (y_[p+1] - y_[p]);
        return (1 - t) * y_[p] + t * y_[p+1] + t * (1-t) * (a * (1-t) + b * t);
    }


    private:
        //find the point to the left
    size_t find_nearest_point( double x ) {
        assert( x >= x_.front() && x <= x_.back() );
        if (x == x_.front())
            return 0;
        if (x == x_.back())
            return x_.size();
        for (size_t r = 1; r < x_.size(); ++r)
            if (x < x_[r])
                return r-1;
        //if we still haven't found anything...
        throw std::out_of_range("couldn't find nearest point... bug in code");
    }



    const std::vector<double>& y_,x_;
    std::vector<double>  coefficient;

};

template <typename scalar>
std::array< scalar, 3 > GCCoefficients( const BasisID &init, const BasisID &fin)
{
    std::array< scalar, 3 > out = {0,0,0};
    //x:
    out[0] = gsl_sf_coupling_3j(init.l*2, 2, fin.l*2, init.m*2, -2, -fin.m*2);
    out[0] -= gsl_sf_coupling_3j(init.l*2, 2, fin.l*2, init.m*2, 2, -fin.m*2);

    out[0] *= gsl_sf_coupling_3j(init.l*2, 2, fin.l*2, 0, 0, 0);
	out[0] *= std::sqrt(( 2 * init.l +1 ) * ( 2 * fin.l +1 ) / 2);
    if ( fin.m%2 == 1 )
        out[0] *= -1;

    //y:
    out[1] = gsl_sf_coupling_3j(init.l*2, 2, fin.l*2, init.m*2, -2, -fin.m*2);
    out[1] += gsl_sf_coupling_3j(init.l*2, 2, fin.l*2, init.m*2, 2, -fin.m*2);

    out[1] *= gsl_sf_coupling_3j(init.l*2, 2, fin.l*2, 0, 0, 0);
	out[1] *= std::complex<scalar>(0, std::sqrt(( 2 * init.l +1 ) * ( 2 * fin.l +1 ) / 2));
    if ( fin.m%2 == 1 )
        out[1] *= -1;

    //z:
    out[2] = gsl_sf_coupling_3j(init.l*2, 2, fin.l*2, 0, 0, 0);
    out[2] *= gsl_sf_coupling_3j(init.l*2, 2, fin.l*2, init.m*2, 0, -fin.m*2);
	out[2] *= std::sqrt(( 2 * init.l +1 ) * ( 2 * fin.l +1 ) );

    return out;
}

void FieldFreePropagate(Vec *H, Vec *wf, PetscReal dt) //propagate forward in time
{
    PetscReal norm1;
    VecNorm(*wf,NORM_2 ,&norm1);

    //copy H:
    Vec tmp;
    VecDuplicate(*H, &tmp);
    VecCopy(*H, tmp);

    //scale by -I:
    VecScale(tmp, std::complex<double>(0,-dt));

    //exponentiate:
    VecExp(tmp);

    //pointwise mult:
    VecPointwiseMult(*wf, *wf, tmp);

    PetscReal norm2;
    VecNorm(*wf,NORM_2 ,&norm2);
    if ( std::abs(norm1 - norm2) / norm1 > 10e-8)
    {
        std::cerr << "before: " << norm1-1 << std::endl;
        std::cerr << "after: " << norm2-1 << " difference: " << norm1 - norm2 <<  std::endl;
    }


    VecDestroy(&tmp);
}


//I should split this into different methods, 
template <typename scalar>
scalar integrateSimpsonsRule(const std::vector<scalar> &psi1, const std::vector<scalar> &psi2, const std::vector<scalar> &grid)
{
    //check for correct size!
    if (psi1.size() != psi2.size() || psi1.size() != grid.size() )
        throw (std::exception());

    std::function < scalar (scalar *) > det = [](scalar* g) {
        return g[0] * g[0] * (g[1] - g[2]) 
             - g[1] * ( g[1] * g[1] - g[2] * g[2] ) 
             + g[1] * g[2] * (g[1] - g[2]);
    };


    std::function< scalar (scalar *, scalar*) > aj =
        [det](scalar* f, scalar* g) {
            scalar aj = f[0] * (g[1] - g[2]) 
                      + f[1] * (g[2] - g[0]) 
                      + f[2] * (g[0] - g[1]);
            aj /= det(g);
            return aj;
        };
    std::function< scalar ( scalar* , scalar*) > bj =
        [det](scalar* f, scalar* g) {
            scalar bj = f[0] * (g[2] * g[2] - g[1] * g[1]) 
                      + f[1] * (g[0] * g[0] - g[2] * g[2]) 
                      + f[2] * (g[0] * g[0] - g[1] * g[1]);
            bj /= det(g);
            return bj;
        };
    std::function< scalar ( scalar* , scalar*) > cj =
        [det](scalar* f, scalar* g) {
            scalar cj = f[0] * g[1] * g[2] * (g[1] - g[2]) 
                      + f[1] * g[0] * g[2] * (g[2] - g[0]) 
                      + f[2] * g[1] * g[0] * (g[0] - g[1]);
            cj /= det(g);
            return cj;
        };
    std::function< scalar ( scalar, scalar, scalar, scalar, scalar, scalar ) > simpsons = 
        [aj, bj, cj]( scalar fjm1, scalar fj, scalar fjp1, scalar gjm1, scalar gj, scalar gjp1 ) {
            scalar g[3] = {gjm1, gj, gjp1};
            scalar f[3] = {fjm1, fj, fjp1};
            scalar result = aj(f, g) * (g[2] * g[2] * g[2] - g[0] * g[0] * g[0]) / 3
                          + bj(f, g) * (g[2] * g[2]        - g[0] * g[0]) / 2
                          + cj(f, g) * (g[2]               - g[0]);
            return result;
        };

    //do the first point, including zero:
    scalar result = simpsons(0., psi1[0] * psi2[0] * grid[0], psi1[1] * psi2[1] * grid[1], 0., grid[0], grid[1]);

    for (size_t i = 1; i < grid.size()-2; i+=2)
    {
        result += simpsons(psi1[i] * psi2[i] * grid[i], psi1[i+1] * psi2[i+1] * grid[i+1] , psi1[i+2] * psi2[i+2] * grid[i+2], grid[i], grid[i+1], grid[i+2]);
    }

    return result;
}

template <typename scalar, typename Func >
scalar integrateTrapezoidRule(const std::vector<scalar> &psi1, const std::vector<scalar> &psi2, const std::vector<scalar> &grid, Func f)
{
    if (psi1.size() != psi2.size() || psi1.size() != grid.size() )
    {
        std::cerr << "integrateTrapezoidRule: the sizes are not the same" << std::endl;
        throw (std::length_error( "integrateTrapezoidRule: the sizes are not the same" ));
    }

    static const std::vector<scalar> dg = [&grid] {
        std::vector<scalar> v;
        v.reserve(grid.size());
        v.push_back(grid[0]/2.);
        for( auto i = begin(grid) + 1; i < end(grid); ++i)
            v.push_back( (*i - *(i - 1))/2.);
        return v;
    }();

    scalar intermediate = psi1[0] * psi2[0] * f(grid[0]);
    scalar result = intermediate * dg[0];
    scalar intermediate1 = 0;

    for (size_t i = 1; i < grid.size(); i++)
    {
        intermediate1 = psi1[i] * f(grid[i]) * psi2[i];
        result += (intermediate1 + intermediate) * dg[i];
        intermediate = intermediate1;
    }

    if (result != result || result == std::numeric_limits<scalar>::infinity() )
    {
        std::cerr << "integrateTrapezoidRule: results in NaN or infinity" << std::endl;
        throw(std::out_of_range( "integrateTrapezoidRule: results in NaN" ));
    }

    return result;
}

//optimized for integrating <psi1|r|psi2>:
template <typename scalar>
scalar integrateTrapezoidRule(const std::vector<scalar> &psi1, const std::vector<scalar> &psi2, const std::vector<scalar> &grid)
{
    //we assume that we only have one grid size that we ever integrate, so we create the difference grid:
    static const std::vector<scalar> dg = [&grid] {
        std::cout << "creating the vector" << std::endl;
        std::vector<scalar> v;
        v.reserve(grid.size());
        v.push_back(grid[0]/2.);
        for( auto i = begin(grid) + 1; i < end(grid); ++i)
            v.push_back( (*i - *(i - 1))/2.);
        return v;
    }();

    scalar intermediate = psi1[0] * psi2[0] * grid[0];
    scalar result = intermediate * dg[0];
    scalar intermediate1 = 0;

    for (size_t i = 1; i < grid.size(); i++)
    {
        intermediate1 = psi1[i] * grid[i] * psi2[i];
        result += (intermediate1 + intermediate) * dg[i];
        intermediate = intermediate1;
    }

    return result;
}

template <typename scalar>
scalar integrateGrid(const std::vector<scalar> &psi1, const std::vector<scalar> &psi2, const std::vector<scalar> &grid)
{
    scalar b = integrateTrapezoidRule(psi1, psi2, grid);
    //scalar a = integrateSimpsonsRule(psi1, psi2, grid);
    //std::cerr << std::scientific;
    //std::cerr.precision(20);
    //std::cerr << "trap: " << a << " simpsons: " << b << " error: " << a-b << " rel: " << a - b / a << std::endl;
    return b;
}

template<typename scalar>
scalar normalize(std::vector<scalar> &wf, const std::vector<scalar> &grid)
{
    //std::function< scalar (scalar) > f = [](scalar r){ return 1.;};
    scalar norm = integrateTrapezoidRule(wf, wf, grid, [](scalar r){return 1.;} );
	norm = std::sqrt(norm);
    for ( auto& a : wf )
        a /= norm;
	return norm;
}

std::vector< double > coulomb_wave_function(kBasisID a, const std::vector<double>& grid)
{
    std::vector< double > cv;
    cv.reserve(grid.size());
    std::complex<double> F, dF;

    Coulomb_wave_functions cwf(true, a.l, -1./a.k);
    for (auto r: grid)
    {
        cwf.F_dF( a.k * r, F, dF);
        cv.push_back(F.real());
    }
    return cv;
}
size_t factorial(size_t n)
{
    if (n == 0)
        return 1;
    return n * factorial(n-1);
};

std::vector<double> gsl_coulomb_wave_function(const kBasisID& a, const std::vector<double>& grid)
{
    int Z = 1;
    std::vector<double> cv;
    cv.reserve(grid.size());
    double* val = new double[1];
    //std::cerr << std::endl << "k: " << a.k << ", l:" << a.l << std::endl;

    for (auto& r : grid)
    {
        double exp;
        //we should really calculate more than just the one l, but not right now.
        gsl_sf_coulomb_wave_F_array ( a.l, 0, - 1./a.k, a.k * r, val, &exp );
        if (*val != *val || *val == std::numeric_limits<double>::infinity() )
        {
            std::cerr << std::endl << "k: " << a.k << ", l:" << a.l << ": " << *val << ", " << r << std::endl;
            throw(std::out_of_range( "gsl_coulomb_wave_function: results in NaN" ));
        }
        cv.push_back(*val);
    }
    normalize(cv, grid); //we want to normalize to 1
    delete[] val;
    return cv;
}

Matrix< double > gsl_spherical_harmonics ( const int lmax, const int n_theta )
{
    Matrix<double> m( n_theta, lmax );

    std::vector< double > t(n_theta);
    int i = 0;
    double dtheta = math::PI / n_theta;
    for ( auto a = m[0]; a < m[n_theta]; a+=lmax )
    {
        gsl_sf_legendre_sphPlm_array(lmax, 0, std::cos( i * dtheta ), a);
        i++;
    }

    return m;
}

inline std::complex<double> Gamma_Lanczos (std::complex<double> z)
{

    std::complex<double> x,t;
    double g = 7.;
    std::vector<double> p(9);
    p[0] = 0.99999999999980993;
    p[1] = 676.5203681218851;
    p[2] = -1259.1392167224028;
    p[3] = 771.32342877765313;
    p[4] = -176.61502916214059;
    p[5] = 12.507343278686905;
    p[6] = -0.13857109526572012;
    p[7] = 9.9843695780195716e-6;
    p[8] = 1.5056327351493116e-7;

    if( real(z) < 0.5)  // reflection
        return PI / (std::sin(PI * z) * Gamma_Lanczos(1. - z));
    else
    {
        z -= 1.;
        x = p[0];
        for( int i=1; i<9; i++) 
            x += p[i]/(z+double(i));
        t = z + g + 0.5;

        return std::sqrt(2.*PI) * std::pow(t,z+0.5) * std::exp(-t) * x;
    }

}

std::tuple<PetscReal, int> VecMax(Vec a)
{
    int loc;
    PetscReal m;
    VecMax(a, &loc, &m);
    return std::make_tuple(m, loc);
}

std::tuple<PetscReal, int> VecMin(Vec a)
{
    int loc;
    PetscReal m;
    VecMin(a, &loc, &m);
    return std::make_tuple(m, loc);
}

std::tuple<PetscReal, int> VecAbsMin(Vec a)
{
    int loc;
    PetscReal m;
    Vec tmp;
    VecDuplicate(a, &tmp);
    VecCopy(a, tmp);
    VecAbs(tmp);
    VecMin(tmp, &loc, &m);
    VecDestroy(&tmp);
    return std::make_tuple(m, loc);
}

//template < typename Func >
//std::vector< std::tuple< PetscScalar, int > > VecSort(Vec a, Func comparator)
//{
    ////mergesort
    //std::vector< int > outint();

template< typename Comparator >
std::vector< std::tuple<PetscScalar, int> > VecFirstNSort(Vec a, size_t n, Comparator comp)
{
    std::vector< PetscInt > outint( n , -1);
    std::vector< PetscScalar > outscalar( n, PetscScalar() );
    std::vector< std::tuple< PetscScalar, int> > out; out.reserve(n);
    PetscScalar *array;
    int low, high;

    VecGetOwnershipRange(a,&low,&high);
    VecGetArray(a, &array);

    MPI_Comm comm;
    int rank, size;
    PetscObjectGetComm((PetscObject)a,&comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);


    //std::cout << *(std::max_element(array, array+(high-low), [](PetscScalar a, PetscScalar b){ return a.real() < b.real(); } )) << std::endl;

    for (int i = 0; i != high-low; ++i)
    {
        auto outn = outscalar.begin();
        auto outm = outint.begin();
        assert (outscalar.end() - outscalar.begin() == outint.end() - outint.begin());
        for (; (outn != outscalar.end()) && (outm != outint.end()) ; ++outm, ++outn)
            if (comp(array[i],*outn))
            {
                outscalar.insert( outn, array[i] );
                outint.insert( outm , i+low );
                outscalar.erase( outscalar.end() - 1 );
                outint.erase( outint.end() - 1 );
                break;
            }
    }

    VecRestoreArray(a, &array);


    //std::cout << "[" << rank << "]" << "send " << std::endl;
    //send to rank 0, with a tag indictating if int or scalar
    if (rank != 0)
    {
        MPI_Send(outscalar.data(), n, MPIU_SCALAR, 0, 0, comm);
        MPI_Send(outint.data(), n, MPIU_INT, 0, 1, comm);
    }

    if (rank == 0)
    {
        std::vector< PetscScalar > allscalars( n * size );
        std::vector< int > allints( n * size );

        //std::cout << "get" << std::endl;
        //loop over all senders:
        for (int senders = 1; senders < size; ++senders)
        {
            MPI_Recv(&allscalars[senders * n], n, MPIU_SCALAR, senders, 0, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&allints[senders * n], n, MPIU_INT, senders, 1, comm, MPI_STATUS_IGNORE);
        }

        std::copy(outscalar.begin(), outscalar.end(), allscalars.begin());
        std::copy(outint.begin(), outint.end(), allints.begin());
        //std::cout << "got" << std::endl;

        //for (auto a : allscalars)
            //std::cout << a << ",";
        //std::cout << std::endl;

        //sort new list:
        for (int i = 0; i != n * size; ++i)
        {
            auto outn = outscalar.begin();
            auto outm = outint.begin();
            assert (outscalar.end() - outscalar.begin() == outint.end() - outint.begin());
            for (; (outn != outscalar.end()) && (outm != outint.end()) ; ++outm, ++outn)
                if (comp(allscalars.at(i),*outn))
                {
                    outscalar.insert( outn, allscalars.at(i) );
                    outint.insert( outm , allints.at(i) );
                    outscalar.erase( outscalar.end() - 1 );
                    outint.erase( outint.end() - 1 );
                    break;
                }
        }

    }
    //std::cout << "[" << rank << "]"<< "rebroadcast " << std::endl;
    MPI_Bcast(&outscalar[0], n, MPIU_SCALAR, 0, comm);
    MPI_Bcast(&outint[0], n, MPIU_INT, 0, comm);

    //std::cout << "got!" << std::endl;
    for (int i = 0; i < n ;  ++i)
        out.push_back(std::make_tuple(outscalar.at(i), outint.at(i)));


    return out;
}

//find the second derivative of a vector, maintaining the vectors size:
template <typename T>
std::vector<T> second_difference( const std::vector<T>& in, const T& h)
{
    std::vector<T> out;
    out.reserve( in.size() );

    //initialize with a forward second derivative:
    out.push_back( (in[2] - 2 * in[1] + in[0]) / (h*h) );

    for (size_t i = 1; i < in.size()-1; ++i)
        out.push_back( (in[i-1] - 2 * in[i] + in[i+1]) / (h*h) );

    out.push_back( (in[ in.size() - 3 ] - 2 * in[in.size() - 2] + in[in.size() - 1]) / (h*h) );

    return out;
}

//fourier
std::vector< std::complex<double> > fourier(std::vector< double >&& time_series )
{
    std::vector< std::complex<double> > out( (time_series.size()) );
    fftw_plan p = fftw_plan_dft_r2c_1d( time_series.size(), 
                                        time_series.data(), 
                                        reinterpret_cast< double(*)[2] >(out.data()), 
                                        FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    return out;
}

//windowing functions:

namespace window
{

    template <typename T>
    std::vector< T > hann_window( std::vector< T >& input )
    {
        //std::vector< T > out( ( input.size() ) );
        for (size_t i = 0; i < input.size(); ++i)
            input[i] *= .5 * (1 - std::cos( 2 * PI * i / (input.size() - 1) ) );

        return input;
    }

}


}
