#include "bspline.hpp"
#include <iostream>
#include <vector>

int main( int arc, char** argv )
{
    std::vector<double>           knots( 10 );
    std::vector<double>::iterator it;
    for ( it = knots.begin(); it < knots.end(); it++ ) *it = it - knots.begin();

    bspline::BSpline<2, double> bs( &knots, 1 );

    for ( int i = 0; i < 100; i++ ) std::cout << bs( i * .1 ) << ",";

    std::cout << std::endl;

    for ( int i = 0; i < 100; i++ )
        std::cout << bs.derivative( i * .1001, 1 ) << ",";

    std::cout << std::endl;
    for ( int i = 0; i < 100; i++ )
        std::cout << bs.derivative( i * .1001, 2 ) << ",";

    return 0;
}
