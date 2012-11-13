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

PetscReal PI = std::atan(1.0)*4.0;
PetscReal C = 137.035999;

template <typename scalar>
scalar CGCoefficient(BasisID init, BasisID fin)
{
    scalar out = 0;
    out = gsl_sf_coupling_3j(init.l*2, 2, fin.l*2, 0, 0, 0);
    out *= out;
	out *= std::sqrt(( 2 * init.l +1 ) * ( 2 * fin.l +1 ) );
    return out;
}

//I should split this into different methods, 
template <typename scalar>
scalar integrateGrid(std::vector<scalar> psi1, std::vector<scalar> psi2, std::vector<scalar> grid)
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



inline std::complex<double> Gamma_Lanczos ( std::complex<double> z)
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
