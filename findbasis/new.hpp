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

void wait_for_key()
{
    std::cout << std::endl << "Press ENTER to continue..." << std::endl;

    std::cin.clear();
    std::cin.ignore( std::cin.rdbuf()->in_avail() );
    std::cin.get();
    return;
}

template <typename scalar>
void display_function( Gnuplot g,
                       std::vector<scalar> grid,
                       std::vector<scalar> wf,
                       int point = -1,
                       bool wait = false,
                       std::array<int, 2> range = {0, 20} )
{
    g.reset_all();
    g.set_xrange(range[0], range[1]);
    g.set_style( "linespoints" )
        .plot_xy( common::vector_type_change<scalar, double>( grid ),
                  common::vector_type_change<scalar, double>( wf ),
                  "wf" );
    if ( point > 0 && point < grid.size() )
        g.set_style( "points" ).plot_xy( std::vector<double>{grid[point]},
                                         std::vector<double>{wf[point]},
                                         "turnover" );
    if ( wait ) wait_for_key();
    g.remove_tempfiles();
}

// the basis to return:
template <typename scalar>
struct basis
{
    std::vector<scalar> wf;
    scalar energy;
    // we also want to keep track of error and error estimates
    scalar de;
    size_t iterations;
};

// the grid description;
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

// describes the iteration
template <typename scalar>
struct iteration
{
    size_t iteration = 0;
    scalar energy_upper;  // Upper bound on the energy
    scalar energy_lower;  // Lower bound on the energy
    bool excited = false; // does the classical turnover point exist?
    scalar energy;        // the energy of this iteration
    int nodes;            // how many nodes for this iteration
    int turnover = -1;    // where the turnover is...

    void upper_bound_bisect() {
        energy_upper = energy;
        energy = (energy_upper + energy_lower) / 2.;
    }

    void lower_bound_bisect() {
        energy_lower = energy;
        energy = (energy_upper + energy_lower) / 2.;
    }
};

template <typename scalar>
void normalize( std::vector<scalar>& wf,
                const std::vector<scalar>& rgrid,
                const scalar dx )
{
    scalar norm = 0;
    for ( size_t i = 0; i < wf.size(); i++ )
        norm += wf[i] * wf[i] * rgrid[i] * rgrid[i] * dx;

    norm = std::sqrt( norm );

    for ( auto& a : wf ) a /= norm;
}

template <typename iterator, typename citerator>
int numerov( citerator fstart,
             citerator fend,
             iterator wfstart,
             iterator wfend )
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

        if ( wfstart[i] < 0 && wfstart[i - 1] >= 0 ) nodes++;
        if ( wfstart[i] > 0 && wfstart[i - 1] <= 0 ) nodes++;
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
std::vector<scalar> make_rgrid( xgrid<scalar> current_grid )
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

//iterator range, exclusive
template <typename iterator>
void print_range( iterator a, iterator b, int star_num = -1 )
{
    for (auto i = a; i < b-1; ++i ) {
        if ( i - a == star_num ) std::cout << "*";
        std::cout << *i << ", ";
    }
    std::cout << b-1 << std::endl << std::endl;
}

template <typename scalar>
int make_f( const std::vector<scalar>& rgrid,
            const BasisID& state,
            const scalar energy,
            const std::function<scalar( scalar, BasisID )>& pot,
            std::vector<scalar>& f )
{
    int turnover;
    f.clear();
    for ( auto r : rgrid ) {
        if ( !f.empty() ) {
            scalar flast = f.back();
            f.push_back(
                    dx * dx / 12 *
                    ( std::pow(
                               ( static_cast<scalar>( state.l ) + .5 ),
                               2 ) +
                      2 * std::pow( r, 2 ) *
                      ( pot( r, state ) - energy ) ) );
            if ( f.back() * flast < 0 && f.back() > 0 )
                turnover = f.size();
        } else {
            f.push_back(
                    dx * dx / 12 *
                    ( std::pow(
                               ( static_cast<scalar>( state.l ) + .5 ),
                               2 ) +
                      2 * std::pow( r, 2 ) *
                      ( pot( r, state ) - energy ) ) );
        }
    }
    for ( auto& fi : f ) fi = 1 - fi;
    return turnover;
}

// perform the numerov from both sides, doesn't normalize because that
// needs rgrid
template <typename scalar>
std::array<int, 2> numerov_from_both_sides( const std::vector<scalar>& f,
                                            std::vector<scalar>& wf,
                                            int turnover )
{
    std::array<int, 2> nodes;
    nodes[0] = numerov( f.begin(),
                        f.begin() + turnover,
                        wf.begin(),
                        wf.begin() + turnover );
    scalar temp = wf[turnover - 1];
    nodes[1] = numerov( f.rbegin() + 1,
                        f.rend() - turnover + 1,
                        wf.rbegin() + 1,
                        wf.rend() - turnover + 1 );
    temp /= wf[turnover - 1];

    // match:
    for ( auto i = wf.rbegin(); i < wf.rend() - turnover + 1; ++i )
        *i *= temp;

    return nodes;
}

template <typename scalar>
int numerov_from_one_side( const std::vector<scalar>& f,
                           std::vector<scalar>& wf )
{
    int nodes;
    nodes = numerov( f.begin(), f.end(), wf.begin(), wf.end() );
    return nodes;
}

template <typename scalar>
std::array<scalar, 2> derivatives( const std::vector<scalar> wf, const std::vector<scalar> f, const int turnover )
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

    auto d21 = forward[0] - 2 * forward[1] + forward[2];
    auto d22 = backward[0] - 2 * backward[1] + backward[2];
    //std::cout << " d2forward: " << d21 << ", d2backward: " << d22
              //<< ", d2fraction: " << d21 / d22 - 1 << std::endl;

    auto d11 = forward[0] - forward[2];
    auto d12 = backward[0] - backward[2];
    //std::cout << " d1forward: " << d11 << ", d1backward: " << d12
              //<< ", d1 subtraction " << d12 - d11 << std::endl;


    return std::array<scalar>{( d12 - d11 ),( d21 / d22 - 1 )};
}


template <typename scalar>
bool find_ground_state( const BasisID state,
                        const xgrid<scalar> grid,
                        const std::function<scalar( scalar, BasisID )> pot,
                        iteration<scalar> it,
                        const scalar err )
{
    // finding the ground state is easier than finding other states.

    Gnuplot g("wf");
    // create a new, smaller grid to iterate on initially:
    auto smallgrid = grid;
    smallgrid.points = static_cast<int>( std::exp( grid.xmax ) );
    auto rgrid = make_rgrid( smallgrid );

    //create the wf and f function from rgrid.
    std::vector<scalar> f(smallgrid.points);
    std::vector<scalar> wf( smallgrid.points,
                            std::numeric_limits<scalar>::quiet_NaN() );

    // find energy_upper.  For the ground state this is the potential at
    // the outside of the grid
    it.energy_upper =
        std::pow( state.l + .5, 2 ) / ( 2 * std::pow( rgrid.back(), 2 ) ) +
        pot( rgrid.back(), state );

    it.energy_lower = 0;
    for ( auto r : rgrid )
        it.energy_lower = std::min( it.energy_lower,
                                    std::pow( state.l + .5, 2 ) /
                                            ( 2. * std::pow( r, 2 ) ) +
                                        pot( r, state ) );
    it.energy = (it.energy_upper + it.energy_lower)/2.;

    // first loop on the smaller grid:
    bool converged = false;
    while ( !converged && it.iteration < 1000 )
    {
        it.turnover = make_f(rgrid, state, pot, it.energy, f);


        if ( it.turnover < 0 || it.turnover > it_grid.points - 2 ) {
            it.turnover = smallgrid.points / 2;
            it.excited = true;
        }
        std::cout << "turnover: " << it.turnover << std::endl;

        // the first couple points:
        wf[0] =
            std::pow( rgrid.front(), state.l + 1 ) *
            ( 1. - 2. * rgrid.front() / ( 2. * scalar( state.l ) + 2. ) ) /
            std::sqrt( rgrid.front() );
        wf[1] = std::pow( rgrid[1], state.l + 1 ) *
                ( 1. - 2. * rgrid[1] / ( 2. * scalar( state.l ) + 2. ) ) /
                std::sqrt( rgrid[1] );

        // the last couple points:
        wf.back() = 0;
        wf[wf.size() - 2] = smallgrid.dx();
        wf[wf.size() - 3] = ( 12 - f[f.size() - 2] * 10 ) *
                            wf[wf.size() - 2] / f[wf.size() - 2];


        auto nodes = numerov_from_both_sides( f, wf, it.turnover );
        normalize( wf, rgrid, smallgrid.dx());


        auto e = derivatives( wf, f, turnover);
        int correct_nodes = state.n - state.l - 1;

        if ( nodes[0] > correct_nodes && e[1] < smallgrid.dx() )
            it.upper_bound_bisect();

        display_function(g, rgrid, wf, it.turnover - 1, true);

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

    tmp = {1, 0, 0, 1, 0};
    iteration<scalar> it;
    scalar err = 1e-18;
    find_ground_state( tmp, desired_grid, pot, it, err );
    //for ( int l = 2; l <= params.lmax(); l++ ) {
        //for ( int n = 20; n <= params.nmax(); n++ ) {
            //tmp = {n, l, 0, 1, 0};
            //find_basis( tmp, desired_grid, pot );
        //}
    //}
}

}

