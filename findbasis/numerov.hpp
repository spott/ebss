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

#include <boost/mpi.hpp>
//#include <boost/serialization.hpp>
//#define DEBUG
//#define DEBUGNODES
//#define DEBUGEXCITED
//#define DEBUGBOUND
//#define DEBUGEND
//#define DEBUGFINAL
// other
#if defined( DEBUG )
#include <../include/gnuplot_i.hpp>
#endif


namespace numerov
{

    std::ofstream err_out;
    std::ofstream reduced_out;

void wait_for_key()
{
    std::cerr << std::endl << "Press ENTER to continue..." << std::endl;

    std::cin.clear();
    std::cin.ignore( std::cin.rdbuf()->in_avail() );
    std::cin.get();
    return;
}

#if defined( DEBUG )
template <typename scalar>
void display_function( Gnuplot& g,
                       const std::vector<scalar>& grid,
                       const std::vector<scalar>& wf,
                       int point = -1,
                       bool wait = false,
                       std::array<double, 2> range = {0, 20},
                       std::string style = "lines" )
{
    g.reset_all();
    g.set_xrange( range[0], range[1] );
    g.set_style( style.c_str() )
        .plot_xy( common::vector_type_change<scalar, double>( grid ),
                  common::vector_type_change<scalar, double>( wf ),
                  "wf" );
    if ( point > 0 && point < grid.size() )
        g.set_style( "points" ).plot_xy(
            std::vector<double>{static_cast<double>( grid[point] )},
            std::vector<double>{static_cast<double>( wf[point] )},
            "turnover" );
    if ( wait ) wait_for_key();
    g.remove_tmpfiles();
}
#endif

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

    scalar dx() const { return ( xmax - xmin ) / ( points - 1 ); }
};

template <typename scalar>
std::ostream& operator<<( std::ostream& out,
                          const xgrid<scalar>& b ) // output
{
    out << "xmin: " << b.xmin << " xmax: " << b.xmax
        << " points: " << b.points << " dx: " << b.dx();
    return out;
}

// describes the iteration
template <typename scalar>
struct iteration
{
    size_t it = 0;
    scalar energy_upper;  // Upper bound on the energy
    scalar energy_lower;  // Lower bound on the energy
    bool excited = false; // does the classical turnover point exist?
    scalar energy;        // the energy of this iteration
    int nodes;            // how many nodes for this iteration
    int turnover = -1;    // where the turnover is...

    void upper_bound_bisect()
    {
        energy_upper = energy;
        energy = ( energy_upper + energy_lower ) / 2.;
    }

    void lower_bound_bisect()
    {
        energy_lower = energy;
        energy = ( energy_upper + energy_lower ) / 2.;
    }
};

template <typename scalar>
std::ostream& operator<<( std::ostream& out,
                          const iteration<scalar>& b ) // output
{
    out << "iteration: " << b.it << " energy_upper: " << b.energy_upper
        << " energy_lower: " << b.energy_lower
        << " de: " << b.energy_upper - b.energy_lower
        << "\n\t energy: " << b.energy << " excited: " << b.excited
        << " turnover: " << b.turnover;
    return out;
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

// iterator range, exclusive
template <typename iterator>
void print_range( iterator a, iterator b, int star_num = -1 )
{
    for ( auto i = a; i < b - 1; ++i ) {
        if ( i - a == star_num ) err_out << "*";
        err_out << *i << ", ";
    }
    err_out << *( b - 1 ) << std::endl << std::endl;
}

template <typename scalar>
int make_f( const std::vector<scalar>& rgrid,
            const BasisID& state,
            const scalar energy,
            const scalar dx,
            const std::function<scalar( scalar, BasisID )>& pot,
            std::vector<scalar>& f )
{
    int turnover = -1;
    f.clear();
    scalar flast = 0;
    for ( const auto r : rgrid ) {
        if ( !f.empty() ) {
            flast = f.back();
            f.push_back(
                dx * dx / 12 *
                ( std::pow( ( static_cast<scalar>( state.l ) + .5 ), 2 ) +
                  2 * std::pow( r, 2 ) * ( pot( r, state ) - energy ) ) );
            if ( f.back() * flast < 0 && f.back() > 0 )
                turnover = f.size();
        } else {
            f.push_back(
                dx * dx / 12 *
                ( std::pow( ( static_cast<scalar>( state.l ) + .5 ), 2 ) +
                  2 * std::pow( r, 2 ) * ( pot( r, state ) - energy ) ) );
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
std::array<scalar, 6> derivatives( const std::vector<scalar> wf,
                                   const std::vector<scalar> f,
                                   const int turnover,
                                   bool errcheck = false )
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

    if ( errcheck ) {
        err_out << "forward: ";
        for ( auto a : forward ) err_out << a << ", ";
        err_out << "\t\tbackward: ";
        for ( auto a : backward ) err_out << a << ", ";
        err_out << std::endl << "wavefunction: ";
        for ( auto a = wf.begin() + turnover - 2;
              a < wf.begin() + turnover + 1;
              ++a )
            err_out << *a << ", ";
        err_out << std::endl;
    }

    auto d21 = forward[0] - 2 * forward[1] + forward[2];
    auto d22 = backward[0] - 2 * backward[1] + backward[2];
    if ( errcheck ) {
        err_out << " d2forward: " << d21 << ", d2backward: " << d22
                  << ", d2fraction: " << d21 / d22 - 1 << std::endl;
    }

    auto d11 = forward[0] - forward[2];
    auto d12 = backward[0] - backward[2];
    if ( errcheck ) {
        err_out << " d1forward: " << d11 << ", d1backward: " << d12
                  << ", d1 subtraction " << d12 - d11 << std::endl;
    }

    auto sign = d11 * d12 < 0 ? math::signum(d11 * forward[0]) : 1;

    return std::array<scalar, 6>{
        { -( d11/d12 - 1 ), ( d21 / d22 - 1 ), d11, d12, d21, d22}};
    //( d21 / d22 - 1 ) * std::abs( d21 - d22 )}};
}


template <typename scalar>
bool
find_ground_state( const BasisID state,
                   const xgrid<scalar> grid,
                   const std::function<scalar( scalar, BasisID )>& pot,
                   iteration<scalar>& it,
                   const scalar err )
{
// finding the ground state is easier than finding other states.
#ifdef DEBUGGROUND
    Gnuplot g( "wf" );
#endif
    err_out << "FINDING GROUND STATE" << std::endl;
    err_out << "====================" << std::endl;
    err_out << "====================" << std::endl;
    // Gnuplot g2( "f" );
    // create a new, smaller grid to iterate on initially:
    auto smallgrid( grid );
    smallgrid.points = static_cast<int>( std::exp( grid.xmax ) );
    auto rgrid = make_rgrid( smallgrid );

    // create the wf and f function from rgrid.
    std::vector<scalar> f( smallgrid.points );
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
    it.energy = ( it.energy_upper + it.energy_lower ) / 2.;

    // first loop on the smaller grid:
    bool converged = false;
    while ( !converged && it.it < 100 ) {
        it.it++;
        it.turnover =
            make_f( rgrid, state, it.energy, smallgrid.dx(), pot, f );
        err_out << "turnover: " << it.turnover << std::endl;
        // display_function( g2,
        // rgrid,
        // f,
        // it.turnover,
        // true,
        //{0, static_cast<double>( rgrid.back() )} );

        if ( it.turnover < 0 || it.turnover > smallgrid.points - 4 ) {
            it.turnover = smallgrid.points / 2;
            it.excited = true;
        }
        err_out << "turnover: " << it.turnover << std::endl;

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
        normalize( wf, rgrid, smallgrid.dx() );


        auto e = derivatives( wf, f, it.turnover );
        int correct_nodes = state.n - state.l - 1;
        it.nodes = nodes[0];

        err_out << " 1st deriv: " << e[0] << ": " << e[2] << "," << e[3] << std::endl
                << " 2nd deriv: " << e[1] << ": " << e[4] << "," << e[5] << std::endl;
        err_out << " nodes: " << nodes[0]
                  << ", correct_nodes: " << correct_nodes << std::endl;
        err_out << it << std::endl;

#ifdef DEBUGGROUND
        display_function(
            g,
            rgrid,
            wf,
            it.turnover - 1,
            true,
            {0., static_cast<double>( rgrid[it.turnover] * 1.5 )} );
#endif

        if ( correct_nodes < it.nodes ) {
            it.upper_bound_bisect();
            continue;
        }
        if ( e[1] > 0 ) {
            it.lower_bound_bisect();
            if ( std::abs( e[1] * (e[2] - e[3]) ) < err ) break;
        }
        if ( e[1] < 0 ) {
            it.upper_bound_bisect();
            if ( std::abs( e[1] * (e[2] - e[3]) ) < err ) break;
        }
    }

    rgrid = make_rgrid( grid );

    // create the wf and f function from rgrid.
    f.resize( grid.points );
    wf = std::vector<scalar>( grid.points,
                              std::numeric_limits<scalar>::quiet_NaN() );

    err_out << "FINDING GROUND STATE pt 2" << std::endl;
    err_out << "====================" << std::endl;
    err_out << "====================" << std::endl;
    while ( !converged && it.it < 1000 ) {
        it.it++;
        it.turnover = make_f( rgrid, state, it.energy, grid.dx(), pot, f );

        if ( it.turnover < 0 || it.turnover > grid.points - 2 ) {
            it.turnover = grid.points / 2;
            it.excited = true;
        }
        err_out << "turnover: " << it.turnover << std::endl;

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
        wf[wf.size() - 2] = grid.dx();
        wf[wf.size() - 3] = ( 12 - f[f.size() - 2] * 10 ) *
                            wf[wf.size() - 2] / f[wf.size() - 2];


        auto nodes = numerov_from_both_sides( f, wf, it.turnover );
        normalize( wf, rgrid, grid.dx() );


        auto e = derivatives( wf, f, it.turnover );
        int correct_nodes = state.n - state.l - 1;
        it.nodes = nodes[0];

        err_out << " 1st deriv: " << e[0] << ": " << e[2] << "," << e[3] << std::endl
                << " 2nd deriv: " << e[1] << ": " << e[4] << "," << e[5] << std::endl;
        err_out << " nodes: " << nodes[0]
                  << ", correct_nodes: " << correct_nodes << std::endl;
        err_out << it << std::endl;

#ifdef DEBUGGROUND
        display_function(
            g,
            rgrid,
            wf,
            it.turnover - 1,
            true,
            {static_cast<double>( rgrid[it.turnover] - 2 ),
             static_cast<double>( rgrid[it.turnover] + 2 )} );
// display_function( g, rgrid, wf, it.turnover - 1, true );
#endif

        if ( it.energy_upper - it.energy_lower < err ) {
            converged = true;
            break;
        }
        if ( correct_nodes < it.nodes ) {
            it.upper_bound_bisect();
            continue;
        }
        scalar d = std::abs(e[1]) > std::abs(e[0]) ? e[1] : e[0];

        if ( d > 0 ) it.lower_bound_bisect();
        if ( d < 0 ) it.upper_bound_bisect();
        //if ( e[1] > 0 ) it.lower_bound_bisect();
        //if ( e[1] < 0 ) it.upper_bound_bisect();
        //if ( e[1] == 0 && e[0] > 0 ) it.lower_bound_bisect();
        //if ( e[1] == 0 && e[0] < 0 ) it.upper_bound_bisect();
        //if ( e[1] == 0 && e[0] == 0 && e[2] - e[3] > 0 ) it.lower_bound_bisect();
        //if ( e[1] == 0 && e[0] == 0 && e[2] - e[3] < 0) it.upper_bound_bisect();
    }
    return converged;
}

template <typename scalar>
bool
converge_on_nodes( const BasisID state,
                   const xgrid<scalar> grid,
                   const std::function<scalar( scalar, BasisID )>& pot,
                   iteration<scalar>& it,
                   const scalar err )
{
// we assume that the energy is at least bounded by the ground state
// and the infinite square well energy

// look up the spherical box solutions for energies... the energies go,
// but there might be some asymptotic solution.
#ifdef DEBUGNODES
    Gnuplot g( "wf" );
    Gnuplot g2( "f" );
#endif

    err_out << "CONVERGE ON NODES" << std::endl;
    err_out << "====================" << std::endl;
    err_out << "====================" << std::endl;
    auto rgrid = make_rgrid( grid );

    // create the wf and f function from rgrid.
    std::vector<scalar> f( grid.points );
    std::vector<scalar> wf( grid.points,
                            std::numeric_limits<scalar>::quiet_NaN() );

    it.energy = ( it.energy_upper + it.energy_lower ) / 2.;
    bool converged = false;

    while ( !converged && it.it < 100 ) {
        it.it++;
        it.turnover = make_f( rgrid, state, it.energy, grid.dx(), pot, f );
        //err_out << "turnover: " << it.turnover << std::endl;
        // display_function( g2,
        // rgrid,
        // f,
        // it.turnover,
        // true,
        //{0, static_cast<double>( rgrid.back() )} );

        if ( it.turnover < 0 || it.turnover > grid.points - 4 ) {
            it.turnover = grid.points / 2;
            it.excited = true;
        } else if ( it.turnover >
                    grid.points * .95 ) { // this is close enough to the
                                          // edge tha twe should be better
                                          // about convergence
            it.excited = true;
        } else
            it.excited = false;
        err_out << "turnover: " << it.turnover << std::endl;

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
        wf[wf.size() - 2] = grid.dx();
        wf[wf.size() - 3] = ( 12 - f[f.size() - 2] * 10 ) *
                            wf[wf.size() - 2] / f[wf.size() - 2];


        if ( !it.excited ) {
            auto n = numerov_from_both_sides( f, wf, it.turnover );
            it.nodes = n[0];
        } else
            it.nodes = numerov_from_one_side( f, wf );
        normalize( wf, rgrid, grid.dx() );

        auto correct_nodes = state.n - state.l - 1;

        err_out << " nodes: " << it.nodes
                  << ", correct_nodes: " << correct_nodes << std::endl;
        err_out << it << std::endl;

#ifdef DEBUGNODES
        display_function( g,
                          rgrid,
                          wf,
                          it.turnover - 1,
                          true,
                          {0., static_cast<double>( rgrid.back() )} );
#endif

        if ( it.nodes < correct_nodes ) {
            it.lower_bound_bisect();
            continue;
        }
        if ( it.nodes > correct_nodes ) {
            it.upper_bound_bisect();
            continue;
        }
        if ( it.nodes == correct_nodes ) {
            if ( it.excited &&
                 ( wf.back() - wf[wf.size() - 2] ) * wf.back() > 0 )
                it.lower_bound_bisect();
            else
                converged = true;
        }
    }
    return converged;
}

template <typename scalar>
bool converge_excited( const BasisID state,
                       const xgrid<scalar> grid,
                       const std::function<scalar( scalar, BasisID )>& pot,
                       iteration<scalar>& it,
                       const scalar err )
{
// we assume that the energy is at least bounded by the ground state
// and the infinite square well energy

// look up the spherical box solutions for energies... the energies go,
// but there might be some asymptotic solution.
#ifdef DEBUGEXCITED
    Gnuplot g( "wf" );
    Gnuplot g2( "f" );
#endif

    err_out << "CONVERGE EXCITED" << std::endl;
    err_out << "====================" << std::endl;
    err_out << "====================" << std::endl;

    auto rgrid = make_rgrid( grid );

    // create the wf and f function from rgrid.
    std::vector<scalar> f( grid.points );
    std::vector<scalar> wf( grid.points,
                            std::numeric_limits<scalar>::quiet_NaN() );

    it.energy = ( it.energy_upper + it.energy_lower ) / 2.;
    bool converged = false;

    int last_nodes = it.nodes;

    while ( !converged && it.it < 1000 ) {
        it.it++;
        it.turnover = make_f( rgrid, state, it.energy, grid.dx(), pot, f );
        err_out << "turnover: " << it.turnover << std::endl;
        // display_function( g2,
        // rgrid,
        // f,
        // it.turnover,
        // true,
        //{0, static_cast<double>( rgrid.back() )} );

        // we already know how this is going to turn out though
        if ( it.turnover < 0 || it.turnover > grid.points - 4 ) {
            it.turnover = grid.points * 3 / 4;
            it.excited = true;
        } else if ( it.turnover < grid.points * .95 ) {
            // if the state is too close to the edge, we want to ignore the
            // "excited" change.
            it.excited = false;
            break;
        } else {
            it.excited = true;
        }
        err_out << "turnover: " << it.turnover << std::endl;

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
        wf[wf.size() - 2] = grid.dx();
        wf[wf.size() - 3] = ( 12 - f[f.size() - 2] * 10 ) *
                            wf[wf.size() - 2] / f[wf.size() - 2];


        if ( !it.excited ) {
            auto n = numerov_from_both_sides( f, wf, it.turnover );
            it.nodes = n[0];
        } else
            it.nodes = numerov_from_one_side( f, wf );
        normalize( wf, rgrid, grid.dx() );

        auto correct_nodes = state.n - state.l - 1;

        err_out << " nodes: " << it.nodes
                  << ", correct_nodes: " << correct_nodes << std::endl;
        err_out << it << std::endl;

#ifdef DEBUGEXCITED
        display_function( g,
                          rgrid,
                          wf,
                          it.turnover - 1,
                          true,
                          {0., static_cast<double>( rgrid.back() )} );
#endif

        if ( it.nodes < correct_nodes ) {
            it.lower_bound_bisect();
            continue;
        }
        if ( it.nodes > correct_nodes ) {
            it.upper_bound_bisect();
            continue;
        }
        if ( it.nodes == correct_nodes ) {
            if ( last_nodes < it.nodes && it.excited )
                it.upper_bound_bisect();
            else if ( last_nodes > it.nodes && it.excited )
                it.lower_bound_bisect();
            else if ( last_nodes == it.nodes )
                it.lower_bound_bisect();
            if ( it.energy_upper - it.energy_lower < err )
                converged = true;
        }
    }
    return converged;
}

template <typename scalar>
bool converge_bound( const BasisID state,
                     const xgrid<scalar> grid,
                     const std::function<scalar( scalar, BasisID )> pot,
                     iteration<scalar>& it,
                     const scalar err )
{
// we assume that the energy is at least bounded by the ground state
// and the infinite square well energy

// look up the spherical box solutions for energies... the energies go,
// but there might be some asymptotic solution.
#ifdef DEBUGBOUND
    Gnuplot g( "wf" );
    Gnuplot g2( "f" );
#endif

    err_out << "CONVERGE BOUND" << std::endl;
    err_out << "====================" << std::endl;
    err_out << "====================" << std::endl;

    auto rgrid = make_rgrid( grid );

    // create the wf and f function from rgrid.
    std::vector<scalar> f( grid.points );
    std::vector<scalar> wf( grid.points,
                            std::numeric_limits<scalar>::quiet_NaN() );

    it.energy = ( it.energy_upper + it.energy_lower ) / 2.;
    bool converged = false;

    while ( !converged && it.it < 1000 ) {
        it.it++;
        it.turnover = make_f( rgrid, state, it.energy, grid.dx(), pot, f );
        //err_out << "turnover: " << it.turnover << std::endl;
        // display_function( g2,
        // rgrid,
        // f,
        // it.turnover,
        // true,
        //{0, static_cast<double>( rgrid.back() )} );

        // we already know how this is going to turn out though
        if ( it.turnover < 0 || it.turnover > grid.points - 4 ) {
            it.turnover = 3 * grid.points / 4;
            it.excited = true;
            break;
        } else
            it.excited = false;
        err_out << "turnover: " << it.turnover << std::endl;

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
        wf[wf.size() - 2] = grid.dx();
        wf[wf.size() - 3] = ( 12 - f[f.size() - 2] * 10 ) *
                            wf[wf.size() - 2] / f[wf.size() - 2];


        if ( !it.excited ) {
            auto n = numerov_from_both_sides( f, wf, it.turnover );
            it.nodes = n[0];
        } else
            it.nodes = numerov_from_one_side( f, wf );
        normalize( wf, rgrid, grid.dx() );

        auto correct_nodes = state.n - state.l - 1;
        auto e = derivatives( wf, f, it.turnover );


        err_out << " 1st deriv: " << e[0] << ": " << e[2] << "," << e[3] << std::endl
                << " 2nd deriv: " << e[1] << ": " << e[4] << "," << e[5] << std::endl;
        err_out << " nodes: " << it.nodes
                  << ", correct_nodes: " << correct_nodes;
        err_out << it << std::endl;

#ifdef DEBUGBOUND
        display_function( g,
                          rgrid,
                          wf,
                          it.turnover - 1,
                          true,
                          {0., static_cast<double>( rgrid.back() )} );
#endif

        if ( it.nodes < correct_nodes ) {
            it.lower_bound_bisect();
            continue;
        }
        if ( it.nodes > correct_nodes ) {
            it.upper_bound_bisect();
            continue;
        }
        if ( it.nodes == correct_nodes ) {
            scalar& d = std::abs(e[1]) > std::abs(e[0]) ? e[1] : e[0];

            if ( d > 0 ) it.lower_bound_bisect();
            if ( d < 0 ) it.upper_bound_bisect();
            //if ( e[1] == 0 && e[0] > 0 ) it.lower_bound_bisect();
            //if ( e[1] == 0 && e[0] < 0 ) it.upper_bound_bisect();
            //if ( e[1] == 0 && e[0] == 0 && e[2] - e[3] > 0 ) it.lower_bound_bisect();
            //if ( e[1] == 0 && e[0] == 0 && e[2] - e[3] < 0) it.upper_bound_bisect();
            if ( it.energy_upper - it.energy_lower < err )
                converged = true;
        }
    }
    return converged;
}

template <typename scalar>
basis<scalar>
find_basis( const BasisID state,
            const xgrid<scalar> grid,
            const std::function<scalar( scalar, BasisID )>& pot,
            scalar err )
{
    bool check = false;
    std::cout << "====================" << std::endl;
    std::cout << "state: " << state << std::endl;
    reduced_out << "====================" << std::endl;
    reduced_out << "state: " << state << std::endl;
    err_out << "====================" << std::endl;
    err_out << "state: " << state << std::endl;

    iteration<scalar> it;
    scalar gserr = 1e-5;
    if ( state.n == 1 ) gserr = err;

    if ( state.e.real() < 0 )
        it.energy_lower = state.e.real();
    else if ( !find_ground_state(
                  {1, 0, 0, 1, 0.}, grid, pot, it, gserr ) ) {
        reduced_out << " couldn't converge ground state! " << std::endl;
        check = true;
    }

    if ( state.n == 1 ) goto converged;


    it.energy_upper =
        std::pow( state.n * math::PI / std::exp( grid.xmax ), 2 ) / 2.;

    if ( !converge_on_nodes(
             state, grid, pot, it, static_cast<scalar>( 1e-8 ) ) ) {
        reduced_out << " couldn't find the node limits! " << std::endl;
    }
    it.it = 0;

    // if not excited before conversion attempt, and not excited after
    // conversion attempt, we failed, otherwise converge as excted state
    if ( !it.excited && !converge_bound( state, grid, pot, it, err ) ) {
        reduced_out << " couldn't converge bound state! " << std::endl;
    }
    if ( it.excited ) {
        it.it = 0;
        if ( !converge_excited( state, grid, pot, it, err ) ) {
            reduced_out << " couldn't converge the excited state! "
                      << std::endl;
            check = it.excited;
        }
    }
    if ( !it.excited ) {
        it.it = 0;
        if ( !converge_bound( state, grid, pot, it, err ) ) {
            reduced_out << " couldn't converge bound state the second time! "
                      << std::endl;
            check = true;
        }
    }
converged:

    err_out << " finished finding energy " << std::endl;
    err_out << it << std::endl;
    auto rgrid = make_rgrid( grid );
    std::vector<scalar> wf( grid.points,
                            std::numeric_limits<scalar>::quiet_NaN() );
    std::vector<scalar> f( grid.points,
                           std::numeric_limits<scalar>::quiet_NaN() );
    auto turnover =
        make_f<scalar>( rgrid, state, it.energy, grid.dx(), pot, f );
#ifdef DEBUG
    Gnuplot g( "final wf" );
#endif

    wf[0] = std::pow( rgrid.front(), state.l + 1 ) *
            ( 1. - 2. * rgrid.front() / ( 2. * scalar( state.l ) + 2. ) ) /
            std::sqrt( rgrid.front() );
    wf[1] = std::pow( rgrid[1], state.l + 1 ) *
            ( 1. - 2. * rgrid[1] / ( 2. * scalar( state.l ) + 2. ) ) /
            std::sqrt( rgrid[1] );

    // the last couple points:
    wf.back() = 0;
    wf[wf.size() - 2] = grid.dx();
    wf[wf.size() - 3] = ( 12 - f[f.size() - 2] * 10 ) * wf[wf.size() - 2] /
                        f[wf.size() - 2];
    bool converged = false;
    while ( !converged && it.it < 1100 ) {
        err_out << it << std::endl;
        err_out << "numerov: " << std::endl;
        numerov_from_both_sides(
            f,
            wf,
            ( std::max(
                it.turnover,
                ( ( turnover > f.size() - 2 ) ? 0 : turnover ) ) ) );
        normalize( wf, rgrid, grid.dx() );
        auto e = derivatives( wf, f, it.turnover );
        err_out << " 1st deriv: " << e[0] << ": " << e[2] << "," << e[3] << std::endl
                << " 2nd deriv: " << e[1] << ": " << e[4] << "," << e[5] << std::endl;
        if (  math::signum( e[2] ) != math::signum( e[3] ) ||
             std::abs( e[0] ) > 0.01 || std::abs( e[1] ) > 0.01 ) {
#ifdef DEBUGEND
            std::cout << " 1st deriv: " << e[0] << ": " << e[2] << "," << e[3] << std::endl
                << " 2nd deriv: " << e[1] << ": " << e[4] << "," << e[5] << std::endl;
            display_function(
                g,
                rgrid,
                wf,
                it.turnover - 1,
                true,
                {static_cast<double>( rgrid[it.turnover] - 2 ),
                 static_cast<double>( rgrid[it.turnover] + 2 )} );
            check=true;
#endif
            it.turnover = ( it.turnover + f.size() - 2 ) / 2;
            it.it++;
        } else {
            converged = true;
        }
    }

    auto e = derivatives( wf, f, it.turnover, false );
    if ( check || std::abs( e[0] ) > 1e-7 || std::abs( e[1] ) > 1e-7 ||
         math::signum( e[2] ) != math::signum( e[3] ) ||
         math::signum( e[4] ) != math::signum( e[5] ) ) {
        e = derivatives( wf, f, it.turnover, true );
        reduced_out << " 1st deriv: " << e[0] << ": " << e[2] << "," << e[3] << std::endl
                << " 2nd deriv: " << e[1] << ": " << e[4] << "," << e[5] << std::endl;
        reduced_out << it << std::endl;
#ifdef DEBUGFINAL
        Gnuplot g( "final wf" );
        Gnuplot g2( "zoomed" );
        display_function( g,
                          rgrid,
                          wf,
                          it.turnover - 1,
                          true,
                          {0, static_cast<double>( rgrid.back() )} );
        display_function(
            g2,
            rgrid,
            wf,
            it.turnover - 1,
            true,
            {static_cast<double>( rgrid[it.turnover] - 1 ),
             static_cast<double>( rgrid[it.turnover] + 1 )} );
#endif
    }

    for ( size_t i = 0; i < wf.size(); i++ ) {
        wf[i] *= std::sqrt( rgrid[i] );
    }
    return basis<scalar>{wf, it.energy, it.energy_upper - it.energy_lower,
                         it.it};
}

// template< typename scalar, typename write_type>
// BasisID basis( BasisID state, const xgrid<scalar> grid, const
// std::function<scalar( scalar, BasisID) > pot, const scalar err)
//{
// auto ret = find_basis( state, grid, pot, err);

// state.e = ret.energy;

// common::export_vector_binary(
// params.basis_function_filename( state ),
// common::vector_type_change<scalar, write_type>( ret.wf ) );

// return state;
//}


template <typename scalar, typename write_type>
void find_basis_set( std::function<scalar( scalar, BasisID )> pot,
                     BasisParameters<scalar, write_type>& params,
                     sae<scalar> atom )
{
#if defined( DEBUG )
    Gnuplot::set_terminal_std( "qt" );
#endif
    xgrid<scalar> desired_grid( {std::log( params.rmin() ),
                                 std::log( params.rmax() ),
                                 static_cast<size_t>( params.points() )} );

    auto rgrid = make_rgrid( desired_grid );
    std::vector<BasisID> energies;
    energies.resize( 0 );
    // params.save_parameters();

    boost::mpi::environment env;
    boost::mpi::communicator world( PETSC_COMM_WORLD,
                                    boost::mpi::comm_attach );

    std::string fname(params.basis_folder() + "/err_log_" +
                 std::to_string( world.rank() ) + ".txt");
    err_out.open(fname , std::ios_base::out);
    fname = std::string(params.basis_folder() + "/reduced_log_" +
                 std::to_string( world.rank() ) + ".txt");
    reduced_out.open(fname , std::ios_base::out);

    basis<scalar> res;
    BasisID tmp;

    energies.resize( 0 );

    iteration<scalar> it;
    scalar err = 1e-15;

    err_out << "Number of threads: " << world.size() << std::endl;

    // find ground state:
    //
    tmp = {1, 0, 0, 0, 0};

    auto ret = find_basis( tmp, desired_grid, pot, err );
    auto gsenergy = ret.energy;

    err_out << "finished finding ground state! " << std::endl;
    // tmp.e = gsenergy;
    // energies.push_back(tmp);
    // std::vector< std::vector< write_type > > output_arrays( std::max(
    // params.lmax() + 1, params.nmax()) );
    if (params.l_only() && params.nmax() - 1 < params.l() )
    {
        std::cerr << " l and nmax don't agree! " << std::endl;
        throw (std::out_of_range(" l and nmax don't agree! "));
    }
    for ( int l = params.l_only() ? params.l() : world.rank();
          l <= (params.l_only()
              ? params.l()
              : std::min( params.lmax(), params.nmax() - 1 ));
          l += world.size() ) {
        err_out << "making output array of size: " << ( params.nmax() - l ) * rgrid.size() << std::endl;
        std::vector<scalar> output_array( ( params.nmax() - l ) *
                                         rgrid.size() );
        std::future<void> block;
        for ( int n = params.n_only() ? params.n() : l + 1; n <= (params.n_only() ? params.n() : params.nmax()); n++ ) {
            err_out << "inner loop: " << n << ", " << l << std::endl;
            if ( n != 1 )
                tmp = {n, l, 0, 0, gsenergy};
            else
                tmp = {n, l, 0, 0, 0};
            ret = find_basis( tmp, desired_grid, pot, err );
            tmp.e = ret.energy;
            // first time through: the future isn't valid.  Second time
            // through, the future is valid,
            // so the first term is false, then we wait, then the second
            // term is true.
            if ( !block.valid() || ( block.wait(), block.valid() ) ) {
                block = std::async(
                    std::launch::async,
                    [n, l, &rgrid, &output_array](
                        std::vector<scalar>&& wf ) {
                        math::Grahm_Schmidt_Orthogonalize(
                            wf,
                            Range<typename std::vector<scalar>::iterator>{
                                output_array.begin(),
                                output_array.begin() +
                                    ( n - ( l + 1 ) ) * wf.size()},
                            rgrid );
                        std::copy( wf.begin(),
                                   wf.end(),
                                   output_array.begin() +
                                       ( n - ( l + 1 ) ) * wf.size() );
                        return;
                    },
                    std::move( ret.wf ) );
            }
            energies.push_back( tmp );
        }

        common::export_vector_binary(
            params.l_block_filename( l ),
            common::vector_type_change<scalar, write_type>(
                output_array ) );
    }

    std::vector<std::vector<BasisID>> all_energies;
    energies.shrink_to_fit();

    if ( world.rank() == 0 )
        gather( world, energies, all_energies, 0 );
    else
        gather( world, energies, 0 );

    energies.clear();

    // auto ei = energies.begin();
    for ( auto i = all_energies.begin(); i < all_energies.end(); ++i ) {
        std::copy_if( i->begin(),
                      i->end(),
                      std::back_inserter( energies ),
                      []( BasisID a ) { return a.n != 0; } );
        // ei += i->size();
    }

    std::sort( energies.begin(), energies.end() );
    if ( world.rank() == 0 )
        for ( auto& a : energies ) reduced_out << a << std::endl;
    std::swap( rgrid, params.grid() );
    std::swap( params.basis_prototype(), energies );
    // params.save_parameters();
}
}

