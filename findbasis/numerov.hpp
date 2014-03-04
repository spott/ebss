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

//#define DEBUGFINAL
// other
#if defined(DEBUG) || defined(DEBUGFINAL)
#include <../include/gnuplot_i.hpp>
#endif


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

#if defined(DEBUG) || defined(DEBUGFINAL)
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
        g.set_style( "points" ).plot_xy( std::vector<double>{grid[point]},
                                         std::vector<double>{wf[point]},
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

    scalar dx() const
    {
        return ( xmax - xmin ) / ( points - 1 );
    }
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

// iterator range, exclusive
template <typename iterator>
void print_range( iterator a, iterator b, int star_num = -1 )
{
    for ( auto i = a; i < b - 1; ++i ) {
        if ( i - a == star_num ) std::cout << "*";
        std::cout << *i << ", ";
    }
    std::cout << *( b - 1 ) << std::endl << std::endl;
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
std::array<scalar, 4> derivatives( const std::vector<scalar> wf,
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
        std::cerr << "forward: ";
        for ( auto a : forward ) std::cerr << a << ", ";
        std::cerr << "\t\tbackward: ";
        for ( auto a : backward ) std::cerr << a << ", ";
        std::cerr << std::endl << "wavefunction: ";
        for ( auto a = wf.begin() + turnover - 2;
              a < wf.begin() + turnover + 1;
              ++a )
            std::cerr << *a << ", ";
        std::cerr << std::endl;
    }

    auto d21 = forward[0] - 2 * forward[1] + forward[2];
    auto d22 = backward[0] - 2 * backward[1] + backward[2];
    if ( errcheck ) {
        std::cerr << " d2forward: " << d21 << ", d2backward: " << d22
                  << ", d2fraction: " << d21 / d22 - 1 << std::endl;
    }

    auto d11 = forward[0] - forward[2];
    auto d12 = backward[0] - backward[2];
    if ( errcheck ) {
        std::cerr << " d1forward: " << d11 << ", d1backward: " << d12
                  << ", d1 subtraction " << d12 - d11 << std::endl;
    }


    return std::array<scalar, 4>{{( d11 - d12 ), ( d21 / d22 - 1 ), d11, d21 }};
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
#ifdef DEBUG
    Gnuplot g( "wf" );
#endif
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
        std::cout << "turnover: " << it.turnover << std::endl;
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
        normalize( wf, rgrid, smallgrid.dx() );


        auto e = derivatives( wf, f, it.turnover );
        int correct_nodes = state.n - state.l - 1;
        it.nodes = nodes[0];

        std::cout << " 1st deriv: " << e[0] << ", 2nd deriv: " << e[1]
                  << std::endl;
        std::cout << " nodes: " << nodes[0]
                  << ", correct_nodes: " << correct_nodes << std::endl;
        std::cout << it << std::endl;

#ifdef DEBUG
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
            if ( std::abs( e[1] * e[0] ) < err ) break;
        }
        if ( e[1] < 0 ) {
            it.upper_bound_bisect();
            if ( std::abs( e[1] * e[0] ) < err ) break;
        }
    }

    rgrid = make_rgrid( grid );

    // create the wf and f function from rgrid.
    f.resize( grid.points );
    wf = std::vector<scalar>( grid.points,
                              std::numeric_limits<scalar>::quiet_NaN() );

    while ( !converged && it.it < 1000 ) {
        it.it++;
        it.turnover = make_f( rgrid, state, it.energy, grid.dx(), pot, f );

        if ( it.turnover < 0 || it.turnover > grid.points - 2 ) {
            it.turnover = grid.points / 2;
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
        wf[wf.size() - 2] = grid.dx();
        wf[wf.size() - 3] = ( 12 - f[f.size() - 2] * 10 ) *
                            wf[wf.size() - 2] / f[wf.size() - 2];


        auto nodes = numerov_from_both_sides( f, wf, it.turnover );
        normalize( wf, rgrid, grid.dx() );


        auto e = derivatives( wf, f, it.turnover );
        int correct_nodes = state.n - state.l - 1;
        it.nodes = nodes[0];

        std::cout << " 1st deriv: " << e[0] << ", 2nd deriv: " << e[1]
                  << std::endl;
        std::cout << " nodes: " << nodes[0]
                  << ", correct_nodes: " << correct_nodes << std::endl;
        std::cout << it << std::endl;

#ifdef DEBUG
        display_function( g,
                          rgrid,
                          wf,
                          it.turnover - 1,
                          true,
                          {static_cast<double>( rgrid[it.turnover] - 2 ),
                           static_cast<double>( rgrid[it.turnover] + 2 )});
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
        if ( e[1] > 0 ) it.lower_bound_bisect();
        if ( e[1] < 0 ) it.upper_bound_bisect();
        if ( e[1] == 0 && e[0] > 0 ) it.lower_bound_bisect();
        if ( e[1] == 0 && e[0] < 0 ) it.upper_bound_bisect();
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
#ifdef DEBUG
    Gnuplot g( "wf" );
    Gnuplot g2( "f" );
#endif

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
        std::cout << "turnover: " << it.turnover << std::endl;
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

        std::cout << " nodes: " << it.nodes
                  << ", correct_nodes: " << correct_nodes << std::endl;
        std::cout << it << std::endl;

#ifdef DEBUG
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
#ifdef DEBUG
    Gnuplot g( "wf" );
    Gnuplot g2( "f" );
#endif

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
        std::cout << "turnover: " << it.turnover << std::endl;
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

        std::cout << " nodes: " << it.nodes
                  << ", correct_nodes: " << correct_nodes << std::endl;
        std::cout << it << std::endl;

#ifdef DEBUG
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
#ifdef DEBUG
    Gnuplot g( "wf" );
    Gnuplot g2( "f" );
#endif

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
        std::cout << "turnover: " << it.turnover << std::endl;
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


        std::cout << " 1st deriv: " << e[0] << ", 2nd deriv: " << e[1]
                  << std::endl;
        std::cout << " nodes: " << it.nodes
                  << ", correct_nodes: " << correct_nodes;
        std::cout << " de: " << it.energy_upper - it.energy_lower
                  << std::endl;
        std::cout << it << std::endl;

#ifdef DEBUG
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
            if ( e[1] > 0 ) {
                it.lower_bound_bisect();
            }
            if ( e[1] < 0 ) {
                it.upper_bound_bisect();
            }
            // else {
            // it.lower_bound_bisect();
            //}
            if ( e[1] == 0 && e[0] > 0 ) {
                it.lower_bound_bisect();
            }
            if ( e[1] == 0 && e[0] < 0 ) {
                it.upper_bound_bisect();
            }
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
    std::cerr << "====================" << std::endl;
    std::cerr << "state: " << state << std::endl;

    iteration<scalar> it;
    scalar gserr = 1e-5;
    if ( state.n == 1 ) gserr = err;

    if ( state.e.real() < 0 )
        it.energy_lower = state.e.real();
    else if ( !find_ground_state(
                   {1, 0, 0, 1, 0.}, grid, pot, it, gserr ) ) {
        std::cerr << " couldn't converge ground state! " << std::endl;
        check = true;
    }

    if ( state.n == 1 ) goto converged;


    it.energy_upper =
        std::pow( state.n * math::PI / std::exp( grid.xmax ), 2 ) / 2.;

    if ( !converge_on_nodes(
              state, grid, pot, it, static_cast<scalar>( 1e-8 ) ) ) {
        std::cerr << " couldn't find the node limits! " << std::endl;
    }
    it.it = 0;

    // if not excited before conversion attempt, and not excited after
    // conversion attempt, we failed, otherwise converge as excted state
    if ( !it.excited && !converge_bound( state, grid, pot, it, err ) ) {
        std::cerr << " couldn't converge bound state! " << std::endl;
        // std::cerr << it << std::endl;
        // std::cerr << state << std::endl;
        // checkA = true;
    }
    if ( it.excited ) {
        it.it = 0;
        if ( !converge_excited( state, grid, pot, it, err ) ) {
            std::cerr << " couldn't converge the excited state! "
                      << std::endl;
            // std::cerr << it << std::endl;
            // std::cerr << state << std::endl;
            check = it.excited;
        }
    }
    if ( !it.excited ) {
        it.it = 0;
        if ( !converge_bound( state, grid, pot, it, err ) ) {
            std::cerr << " couldn't converge bound state the second time! "
                      << std::endl;
            // std::cerr << it << std::endl;
            // std::cerr << state << std::endl;
            check = true;
        }
    }
converged:

    std::cout << it << std::endl;
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
#ifdef DEBUG
    display_function( g,
                      rgrid,
                      f,
                      it.turnover - 1,
                      true,
                      {static_cast<double>( rgrid[it.turnover] - 2 ),
                       static_cast<double>( rgrid[it.turnover] + 2 )}
                      );
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
        std::cout << it << std::endl;
        numerov_from_both_sides(
            f,
            wf,
            ( std::max(
                it.turnover,
                ( ( turnover > f.size() - 2 ) ? 0 : turnover ) ) ) );
        normalize( wf, rgrid, grid.dx() );
        auto e = derivatives( wf, f, it.turnover );
        if ( std::abs(e[3] * e[0]) > 0.01  || std::abs( e[1] * e[4] ) > 0.01 ) {
#ifdef DEBUG
            display_function(
                g,
                rgrid,
                wf,
                it.turnover - 1,
                true,
                {static_cast<double>( rgrid[it.turnover] - 2 ),
                 static_cast<double>( rgrid[it.turnover] + 2 )});
#endif

            it.turnover = ( it.turnover + f.size() - 2 ) / 2;
            it.it++;
        } else {
            converged = true;
        }
    }

    auto e = derivatives( wf, f, it.turnover, false );
    if ( true || std::abs( e[0] ) > .1 || std::abs( e[1] ) > .1 ) {
        e = derivatives( wf, f, it.turnover, true );
        std::cerr << " 1st deriv: " << e[0] << " 2nd deriv: " << e[1]
                  << std::endl;
        std::cerr << it << std::endl;
#ifdef DEBUGFINAL
        Gnuplot g( "final wf" );
        Gnuplot g2( "zoomed" );
        display_function( g,
                          rgrid,
                          wf,
                          it.turnover - 1,
                          true,
                          {0, static_cast<double>( rgrid.back() )} );
        display_function( g2,
                          rgrid,
                          wf,
                          it.turnover - 1,
                          true,
                          {static_cast<double>( rgrid[it.turnover] - 1 ),
                           static_cast<double>( rgrid[it.turnover] + 1 )});
#endif
    }

    for ( size_t i = 0; i < wf.size(); i++ ) {
        wf[i] *= std::sqrt( rgrid[i] );
    }
    return basis<scalar>{wf,                                it.energy,
                         it.energy_upper - it.energy_lower, it.it};
}

//template< typename scalar, typename write_type>
//BasisID basis( BasisID state, const xgrid<scalar> grid, const std::function<scalar( scalar, BasisID) > pot, const scalar err)
//{
    //auto ret = find_basis( state, grid, pot, err);

    //state.e = ret.energy;

    //common::export_vector_binary(
            //params.basis_function_filename( state ),
            //common::vector_type_change<scalar, write_type>( ret.wf ) );

    //return state;
//}


template <typename scalar, typename write_type>
void find_basis_set( std::function<scalar( scalar, BasisID )> pot,
                     BasisParameters<scalar, write_type>& params,
                     sae<scalar> atom )
{
#if defined(DEBUG) || defined(DEBUGFINAL)
    Gnuplot::set_terminal_std( "qt" );
#endif
    xgrid<scalar> desired_grid( {std::log( params.rmin() ),
                                 std::log( params.rmax() ),
                                 static_cast<size_t>( params.points() )} );

    auto rgrid = make_rgrid(desired_grid);
    std::swap(rgrid,params.grid());
    std::vector<BasisID> energies;
    energies.resize(0);
    //params.save_parameters();

    basis<scalar> res;
    BasisID tmp;

    energies.resize( 0 );

    iteration<scalar> it;
    scalar err = 1e-15;

    size_t num_threads = params.procs();
    std::cout << "Number of threads: " << num_threads << std::endl;

    // find ground state:
    tmp = {1, 0, 0, 0, 0};
    auto ret = find_basis( tmp, desired_grid, pot, err );
    auto gsenergy = ret.energy;
    //tmp.e = gsenergy;
    //energies.push_back(tmp);
    //std::vector< std::vector< write_type > > output_arrays( std::max( params.lmax() + 1, params.nmax()) );

    for ( int l = 0; l <= std::min(params.lmax(), params.nmax()-1); l++ ) {
        std::vector< scalar > output_array( ( params.nmax() - l) * params.grid().size() );
        for ( int n = l + 1; n <= params.nmax(); n++ ) {
            if (n != 1)
                tmp = {n, l, 0, 0, gsenergy};
            else
                tmp = {n, l, 0, 0, 0};
            ret = find_basis( tmp, desired_grid, pot, err );
            tmp.e = ret.energy;
            std::cerr << " orthogonalizing! " << std::endl;
            math::Grahm_Schmidt_Orthogonalize(
                ret.wf,
                Range<typename std::vector<scalar>::iterator>{
                    output_array.begin(),
                    output_array.begin() +
                        ( n - ( l + 1 ) ) * ret.wf.size()},
                params.grid() );
            std::copy( ret.wf.begin(),
                       ret.wf.end(),
                       output_array.begin() +
                           ( n - ( l + 1 ) ) * ret.wf.size() );
            energies.push_back(tmp);
        }
        common::export_vector_binary( params.l_block_filename(l), common::vector_type_change<scalar, write_type>(output_array) );
    }

    std::sort( energies.begin(), energies.end() );
    for ( auto& a : energies ) std::cout << a << std::endl;
    std::swap( params.basis_prototype(), energies);
    //params.save_parameters();
}
}

