#pragma once
//ebss
#include<common/parameters/Parameters.hpp>

//stl:
#include<complex>
#include<vector>
#include<cmath>

//gsl:
#include<gsl/gsl_sf_coupling.h>

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
    scalar out = 0;
    out = gsl_sf_coupling_3j(init.l*2, 2, fin.l*2, 0, 0, 0);
    out *= out;
	out *= std::sqrt(( 2 * init.l +1 ) * ( 2 * fin.l +1 ) );
    return out;
}

void FieldFreePropagate(Vec *H, Vec *wf, PetscReal dt) //propagate forward in time
{
    PetscReal norm1;
    VecNorm(*wf,NORM_2 ,&norm1);
    std::cout << "before: " << norm1-1 << std::endl;

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
    std::cout << "after: " << norm2-1 << " difference: " << norm1 - norm2 <<  std::endl;


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


template <typename scalar>
scalar integrateTrapezoidRule(const std::vector<scalar> &psi1, const std::vector<scalar> &psi2, const std::vector<scalar> &grid)
{
    //check for correct size!
    if (psi1.size() != psi2.size() || psi1.size() != grid.size() )
        throw (std::exception());

    scalar result = psi1[0] * psi2[0] * grid[0] * grid[0] / 2;

    for (size_t i = 1; i < grid.size(); i++)
        result += (psi1[i] * grid[i] * psi2[i] + psi1[i-1] * grid[i-1] * psi2[i-1]) * (grid[i] - grid[i-1]) / 2;

    if (result != result || result == std::numeric_limits<scalar>::infinity() )
        throw(std::exception());

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
	//auto wf = *wavefunction;
	scalar norm = wf[0] * wf[0] * grid[0] / 2; //the first point is always zero...
	scalar tmp1 = norm;
	scalar tmp2 = 0;
	for (size_t i = 1; i < grid.size(); i++)
	{
		tmp2 = wf[i] * wf[i] * (grid[i]-grid[i-1]) / 2;
		norm += tmp2 + tmp1;
		tmp1 = tmp2;
	}

    if (norm != norm || norm == std::numeric_limits<scalar>::infinity() )
        throw(std::exception());
	norm = std::sqrt(norm);
	for (size_t i = 0; i < grid.size(); i++)
	{
		wf[i] /= norm;
	}
	return norm;
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
}
