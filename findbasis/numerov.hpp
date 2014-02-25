#pragma once

// stl
#include <vector>
#include <forward_list>
#include <array>
#include <iterator>
#include <cmath>
#include <iostream>
#include <limits>
#include <future>
#include <stdexcept>

// ebss
#include <common/common.hpp>
#include <common/math.hpp>
#include <findbasis/single_active_electron.hpp>

// other
#include <../include/gnuplot_i.hpp>

namespace numerov
{

    //forward declare:

template <typename iterator, typename citerator>
int
numerov( citerator fstart, citerator fend, iterator wfstart, iterator wfend );

void wait_for_key()
{
    std::cout << std::endl << "Press ENTER to continue..." << std::endl;

    std::cin.clear();
    std::cin.ignore( std::cin.rdbuf()->in_avail() );
    std::cin.get();
    return;
}

template <typename scalar>
struct basis
{
    std::vector<scalar> wf;
    scalar energy;
    // we also want to keep track of error and error estimates
    scalar de;
    size_t iterations;
};

template <typename scalar>
struct xgrid
{
    scalar xmin, xmax;
    size_t points;

    scalar dx()
    {
        return ( xmax - xmin ) / ( points - 1 );
    }
};

template <typename scalar>
struct it
{
    size_t iteration = 0;
    scalar energy_upper;  // Eventually will be the n-l-1 == 1 energy
    scalar energy_lower;  // This won't change, except for special cases
                          // (but we
                          // will test for it anyways...
    bool excited = false; // does the classical turnover point exist?
    scalar energy;        // the energy of this iteration
    int nodes;            // how many nodes for this iteration
    int turnover = -1;    // where the turnover is...
    scalar f;             // the derivative matching value;
    scalar de;            // Will be used to bisect
                          // bool excited;
};

template <typename scalar>
std::ostream& operator<<( std::ostream& out,
                          const it<scalar>& b ) // output
{
    out << "iteration: " << b.iteration
        << " energy_upper: " << b.energy_upper
        << " energy_lower: " << b.energy_lower
        << "\n\t energy: " << b.energy << " f: " << b.f << " de: " << b.de;
    return out;
}

template <typename scalar>
scalar forward_derivative( const std::vector<scalar>& wf, const std::vector<scalar>& r, size_t p)
{
    return (wf[p] - wf[p-1]) / (r[p] - r[p-1]);
}

template <typename scalar>
scalar derivative2( const std::vector<scalar>& wf, const std::vector<scalar>& r, size_t p)
{
    return (wf[p-1] - 2. * wf[p] + wf[p+1]) / ((r[p] - r[p-1])*(r[p+1] - r[p]));
}

// 1st derivative check:
template <typename scalar>
scalar d1_check( const std::vector<scalar> wf,
                 const std::vector<scalar> f,
                 const int turnover )
{
    // copy the wavefunction, and then propagate it around the turnover
    // point (both forward and backward)
    std::vector<scalar> forward( 3 );
    std::copy( wf.begin() + turnover - 2,
               wf.begin() + turnover + 1,
               forward.begin() );
    numerov( f.begin() + turnover - 2,
             f.begin() + turnover + 1,
             forward.begin(),
             forward.end() );
    std::vector<scalar> backward( 3 );
    std::copy( wf.rend() - turnover - 1,
               wf.rend() - turnover + 2,
               backward.rbegin() );
    numerov( f.rend() - turnover - 1,
             f.rend() - turnover + 2,
             backward.rbegin(),
             backward.rend() );

    //std::cout << " wf " << std::endl;
    //for(auto a = wf.begin() + turnover - 2; a < wf.begin() + turnover + 1; ++a)
        //std::cout << *a << ", ";
    //std::cout << std::endl << " forward: " << std::endl;
    //for(auto a: forward)
        //std::cout << a << ", ";
    //std::cout << std::endl << " backward: " << std::endl;
    //for(auto a: backward)
        //std::cout << a << ", ";
    //std::cout << std::endl;

    auto d21 = forward[0] - 2 * forward[1] + forward[2];
    auto d22 = backward[0] - 2 * backward[1] + backward[2];
    std::cout << " d2forward: " << d21 << ", d2backward: " << d22 << ", d2fraction: " << d21/d22 - 1 << std::endl;

    auto d11 = forward[0] -  forward[2];
    auto d12 = backward[0] - backward[2];
    std::cout << " d1forward: " << d11 << ", d1backward: " << d12 << ", d1 subtraction " << d12 - d11 << std::endl;


    return std::abs((d12 - d11)) * (d21 / d22 - 1);

}

// check for poor wf by checking for "smoothness":
template <typename scalar>
scalar max_deriv( const std::vector<scalar>& wf,
                  const std::vector<scalar>& rgrid )
{
    scalar max_deriv = 0;
    for ( size_t i = 1; i < wf.size(); ++i )
        max_deriv = std::max( std::abs( ( wf[i] - wf[i - 1] ) /
                                        ( rgrid[i] - rgrid[i - 1] ) ),
                              max_deriv );
    return max_deriv;
}

// auxillary functions:

template <typename scalar>
bool ordered( scalar v1, scalar v2, scalar v3)
{
    return (v1 > v2 && v2 > v3) || (v1 < v2 && v2 < v3);
}


template <typename scalar>
void normalize( std::vector<scalar>& wf,
                const std::vector<scalar>& rgrid,
                const scalar dx )
{
    scalar norm = 0;
    for ( size_t i = 0; i < wf.size(); i++ )
        norm += wf[i] * wf[i] * rgrid[i] * rgrid[i] * dx;

    norm = std::sqrt( norm );

    for ( auto& a : wf ) 
        a /= norm;
}


template <typename scalar>
void print_vector( std::vector<scalar> in, int star_num = -1 )
{
    for ( auto i = in.begin(); i < in.end() - 1; ++i ) {
        if ( i - in.begin() == star_num ) std::cout << "*";
        std::cout << *i << ", ";
    }
    std::cout << in.back() << std::endl << std::endl;
}


template <typename scalar>
std::vector<scalar> make_rgrid( scalar dx, xgrid<scalar> current_grid )
{
    std::vector<scalar> g;
    g.reserve( current_grid.points );
    for ( size_t i = 0; i < current_grid.points; ++i )
        g.push_back(
            std::exp( current_grid.xmin +
                      i * ( current_grid.xmax - current_grid.xmin ) /
                          ( current_grid.points - 1 ) ) );
    return g;
}


template <typename iterator, typename citerator>
int
numerov( citerator fstart, citerator fend, iterator wfstart, iterator wfend )
{
    unsigned int nodes = 0;

    if ( wfend - wfstart != fend - fstart || wfend - wfstart == 0 )
        std::cerr << "wf and f don't have the same size or are zero."
                  << std::endl;
    std::cout << std::endl;

    bool inf = false;
    for ( int i = 2; i < wfend - wfstart; i++ ) {
        if ( !inf )
            wfstart[i] = ( ( 12. - fstart[i - 1] * 10. ) * wfstart[i - 1] -
                           fstart[i - 2] * wfstart[i - 2] ) /
                         fstart[i];
        else
            wfstart[i] = wfstart[i - 1];

        //std::cout << wfstart[i] << ", " << wfstart[i-1] << ": " <<  (wfstart[i] < 0) << ", " << (wfstart[i-1] >= 0) << std::endl;
        if ( wfstart[i] < 0 && wfstart[i-1] >= 0 ) {
            nodes++;
            //std::cout << " node detected!" << std::endl;
        }
        if ( wfstart[i] > 0 && wfstart[i-1] <= 0 ) {
            nodes++;
            //std::cout << " node detected!" << std::endl;
        }
        // check for inf & nan:
        if ( std::abs( *( wfstart + i ) ) ==
                 std::numeric_limits<
                     typename iterator::value_type>::infinity() ||
             ( wfstart[i] ) != ( wfstart[i] ) ) {
            wfstart[i] = wfstart[i - 1];
            std::cerr << "infinity detected at " << i
                      << " last two values: " << *( wfstart + i - 1 )
                      << ", " << wfstart[i - 2] << " nodes: " << nodes
                      << std::endl;
            throw( std::out_of_range( "infinity detected" ) );
            inf = true;
        }
    }

    return nodes;
}

template <typename scalar>
basis<scalar>
find_basis( const BasisID state,
            const xgrid<scalar> desired_grid, // desired xgrid
            const std::function<scalar( scalar, BasisID )> pot )
{
    std::cout << "state: " << state << std::endl;
    // start with a small grid, iterate until we get a decent convergence,
    // then embiggen the grid and repeat till we are at
    // our desired grid size.

    // the starting grid should be a power of two time smaller than the
    // desired (so I can just double it as I move forward) but we need to
    // know how small to go.  Since we are on a logarithmic grid, and dx is
    // the difference in the derivative, we will make a choice that dx is <
    // .5 for our smallest grid:
    const scalar dx_target = .1;
    const scalar maximum_derivative = 1000.;

    xgrid<scalar> current_grid = desired_grid;
    while ( ( current_grid.points /= 2 ) >
            ( ( current_grid.xmax - current_grid.xmin ) / dx_target ) ) {
    }

    Gnuplot g1( "wf" );
    Gnuplot g2( "f" );

    //setup the initial stuff:
    const scalar dx = current_grid.dx();
    const std::vector<scalar> rgrid = make_rgrid( dx, current_grid );
    it<scalar> current;


    // we are aloways lower than a particle in a box, so the "energy
    // upper" is definitely upper bound by a particle in a box, where
    // the energy is E = n^2 \pi^2 / 2 L^2 (in atomic units)
    current.energy_upper =
        std::pow( state.n * math::PI / rgrid.back(), 2 ) / 2.;

    // we are also always higher than the bottom of the potential well
    current.energy_lower = 0;
    for ( auto r : rgrid ) {
        current.energy_lower =
            std::min( current.energy_lower,
                      std::pow( scalar( state.l ) + .5, 2 ) /
                              ( 2. * std::pow( r, 2 ) ) +
                          pot( r, state ) );
    }
    // the energy guess is an average of the two:
    current.energy =
        ( current.energy_upper + current.energy_lower ) / 2.;
    //current.energy = -.5;

    // loop until the number of points is correct
    while ( current_grid.points <= desired_grid.points ) {
        std::cout << "current grid points: " << current_grid.points
                  << "current grid dx: " << current_grid.dx()
                  << std::endl;
        // create the rgrid, find the dx:
        const scalar dx = current_grid.dx();
        const std::vector<scalar> rgrid = make_rgrid( dx, current_grid );
        //std::cout << "rgrid: " << std::endl;
        //print_vector( rgrid );
        // the function and wf:
        std::vector<scalar> f;
        f.reserve( current_grid.points );
        std::vector<scalar> wf( current_grid.points,
                                std::numeric_limits<scalar>::quiet_NaN() );




        // we want to be able to traverse our limits in not too many steps.
        current.de = ( current.energy_upper - current.energy_lower ) / 10.;

        bool converged = false;

        while ( !converged && current.iteration < 10000 ) {
            std::cout << current << std::endl;
            current.iteration++;

            // this probably isn't the best way to do this...
            f.clear();
            current.turnover = -1;
            // populate the f vector and figure out where the turnover
            // point is.
            for ( auto r : rgrid ) {
                if ( !f.empty() ) {
                    scalar flast = f.back();
                    f.push_back(
                        dx * dx / 12 *
                        ( std::pow(
                              ( static_cast<scalar>( state.l ) + .5 ),
                              2 ) +
                          2 * std::pow( r, 2 ) *
                              ( pot( r, state ) - current.energy ) ) );
                    if ( f.back() * flast < 0 && f.back() > 0 )
                        current.turnover = f.size();
                } else {
                    f.push_back(
                        dx * dx / 12 *
                        ( std::pow(
                              ( static_cast<scalar>( state.l ) + .5 ),
                              2 ) +
                          2 * std::pow( r, 2 ) *
                              ( pot( r, state ) - current.energy ) ) );
                }
            }
            for ( auto& fi : f ) fi = 1 - fi;
            //std::cout << "f: " << std::endl;
            //print_vector( f, current.turnover );

            // if the turnover is less than zero (meaning we didn't find
            // one, or the turnover isn't enough points from the end, then
            // we want to evenly split the grid and move from both ends at
            // once.
            //
            // TODO we might want to keep track of if we do this check or
            // not.
            if ( current.turnover < 0 ||
                 current.turnover > current_grid.points - 2 ) {
                current.turnover = current_grid.points / 2;
                current.excited = true;
            }
            std::cout << "turnover: " << current.turnover << " excited!: " << current.excited << std::endl;

            // the first couple points:
            wf[0] = std::pow( rgrid.front(), state.l + 1 ) *
                    ( 1. - 2. * rgrid.front() /
                               ( 2. * scalar( state.l ) + 2. ) ) /
                    std::sqrt( rgrid.front() );
            wf[1] =
                std::pow( rgrid[1], state.l + 1 ) *
                ( 1. - 2. * rgrid[1] / ( 2. * scalar( state.l ) + 2. ) ) /
                std::sqrt( rgrid[1] );

            wf.back() = 0;
            wf[wf.size() - 2] = dx;
            wf[wf.size() - 3] = ( 12 - f[f.size() - 2] * 10 ) *
                                wf[wf.size() - 2] / f[wf.size() - 2];


            // all derivatives should be continuous, so we look at
            // current.turnover
            current.nodes = numerov( f.begin(),
                                     f.begin() + current.turnover,
                                     wf.begin(),
                                     wf.begin() + current.turnover );
            auto nodes_before = current.nodes;
            scalar temp = wf[current.turnover - 1];
            std::cout << "midpoint from one side: " << temp;
            std::cout << " nodes_before: " << nodes_before;
            current.nodes += numerov( f.rbegin() + 1,
                                      f.rend() - current.turnover + 1,
                                      wf.rbegin() + 1,
                                      wf.rend() - current.turnover + 1 );
            temp /= wf[current.turnover-1];
            std::cout << " total nodes: " << current.nodes;
            if ( !current.excited ) 
                current.nodes = nodes_before;
            std::cout << " actual nodes: " << current.nodes << std::endl;
            // normalize:
            for ( auto i = wf.rbegin(); i < wf.rend() - current.turnover + 1;
                  ++i )
                *i *= temp;

            normalize( wf, rgrid, dx );

            // check the max derivative:

            //std::cout << "wf: " << std::endl;
            //print_vector( wf, current.turnover - 1 );




            auto m = max_deriv( wf, rgrid );
            std::cout << "max deriv: " << m << std::endl;
            // the checking code:

            //check the derivatives at the cusp:
            g1.reset_all();
            //g1.set_xrange(0,20);
            g1.set_style( "linespoints" ).plot_xy( common::vector_type_change<scalar,double>(rgrid), common::vector_type_change<scalar,double>(wf), "wf" );
            g1.set_style( "points" )
                .plot_xy( std::vector<double>{rgrid[current.turnover-1]},
                          std::vector<double>{wf[current.turnover-1]},
                          "turnover" );
            g2.reset_all();
            //g2.set_xrange(0,20);
            g2.set_style( "linespoints" ).plot_xy( common::vector_type_change<scalar,double>(rgrid), common::vector_type_change<scalar,double>(f), "f" );
            g2.set_style( "points" )
                .plot_xy( std::vector<double>{rgrid[current.turnover-1]},
                          std::vector<double>{f[current.turnover-1]},
                          "turnover" );
            g2.set_style("linespoints" ).plot_xy( std::vector<double>{0,10}, std::vector<double>{1,1});
            wait_for_key();
            g1.remove_tmpfiles();
            g2.remove_tmpfiles();

            //check the nodes:


            auto e = d1_check( wf, f, current.turnover);

            if ( current.nodes > state.n - state.l - 1 &&
                 std::abs( e ) < current_grid.dx() ) {
                current.energy_upper = current.energy;
                current.energy =
                    ( current.energy_lower + current.energy_upper ) / 2;
                continue;
            }
            if ( current.nodes < state.n - state.l - 1 &&
                 std::abs( e ) < current_grid.dx() ) {
                current.energy_lower = current.energy;
                current.energy =
                    ( current.energy_lower + current.energy_upper ) / 2;
                continue;
            }
            std::cout << "e: " << e << " right nodes " << std::endl;

            if ( e > 0 ) {
                current.energy_lower = current.energy;
                current.energy =
                    ( current.energy_lower + current.energy_upper ) / 2;
                if ( std::abs( e ) < current_grid.dx()/ 2. &&
                     current_grid.points * 2 < desired_grid.points )
                    break;
            }
            if ( e < 0 ) {
                if ( current.nodes > state.n - state.l - 1 )
                    current.energy_upper = current.energy;
                current.energy -= current.de;
                if ( std::abs( e ) < current_grid.dx()/2. &&
                     current_grid.points * 2 < desired_grid.points )
                    break;
            }
            std::cout << "===============================" << std::endl;
        }
        current_grid.points *= 2;
    }
}

template <typename scalar, typename write_type>
void find_basis_set( std::function<scalar( scalar, BasisID )> pot,
                     BasisParameters<scalar, write_type>& params,
                     sae<scalar> atom )
{
    Gnuplot::set_terminal_std( "qt" );
    xgrid<scalar> desired_grid( {std::log( params.rmin() ),
                                 std::log( params.rmax() ),
                                 static_cast<size_t>( params.points() )} );

    std::vector<BasisID> energies = params.basis_prototype();
    params.save_parameters();

    basis<scalar> res;
    BasisID tmp;

    energies.resize( 0 );

    for ( int l = 2; l <= params.lmax(); l++ ) {
        for ( int n = 30; n <= params.nmax(); n++ ) {
            tmp = {n, l, 0, 1, 0};
            find_basis( tmp, desired_grid, pot );
        }
    }
}
}
