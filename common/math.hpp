#include<common/parameters.hpp>

namespace math{
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
