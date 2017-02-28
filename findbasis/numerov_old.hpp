#pragma once

// stl
#include <array>
#include <cmath>
#include <forward_list>
#include <future>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <vector>

// ebss
#include <common/common.hpp>
#include <common/math.hpp>
#include <findbasis/single_active_electron.hpp>

namespace numerov
{
template <typename iterator>
int numerov( iterator fstart, iterator fend, iterator wfstart, iterator wfend )
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

        // std::cerr << i << " " << wfstart[i] << " ";
        if ( *( wfstart + i ) < 0 && *( wfstart + i - 1 ) >= 0 ) nodes++;
        if ( *( wfstart + i ) > 0 && *( wfstart + i - 1 ) <= 0 ) nodes++;
        // check for inf & nan:
        if ( std::abs( *( wfstart + i ) ) ==
                 std::numeric_limits<
                     typename iterator::value_type>::infinity() ||
             ( wfstart[i] ) != ( wfstart[i] ) ) {
            wfstart[i] = wfstart[i - 1];
            std::cerr << "infinity detected at " << i
                      << " last two values: " << *( wfstart + i - 1 ) << ", "
                      << wfstart[i - 2] << " nodes: " << nodes << std::endl;
            throw( std::out_of_range( "infinity detected" ) );
            inf = true;
        }
    }

    return nodes;
}

template <typename scalar>
struct basis {
    std::vector<scalar> wf;
    scalar              energy;
};
template <typename scalar>
struct it {
    int    iteration;
    scalar energy_upper;   // Eventually will be the n-l-1 == 1 energy
    bool upper_converged;  // this tells us if the energy_upper is converged on
                           // the correct value;
    scalar energy_lower;  // This won't change, except for special cases (but we
                          // will test for it anyways...
    scalar energy;        // the energy of this iteration
    int    nodes;         // how many nodes for this iteration
    int    turnover;      // where the turnover is...
    scalar f;             // the derivative matching value;
    scalar de;            // Will be used to bisect
    // bool excited;
};


template <typename scalar>
std::ostream& operator<<( std::ostream& out, const it<scalar>& b )  // output
{
    out << "iteration: " << b.iteration << " energy_upper: " << b.energy_upper
        << " upper_converged: " << ( b.upper_converged ? "true" : "false" )
        << " energy_lower: " << b.energy_lower << "\n\t energy: " << b.energy
        << " f: " << b.f << " de: " << b.de;
    return out;
};

template <typename scalar>
basis<scalar> find_basis( const BasisID state, const scalar dx,
                          const std::vector<scalar>& rgrid,
                          const std::function<scalar( scalar, BasisID )> pot,
                          const scalar err )
{
    int messiness =
        3000;  // where the messiness for the shitty j's is adjusted away:
    std::array<scalar, 3> fore;  // I don't think I need a memory of these...
    std::array<scalar, 3> back;
    basis<scalar> result;  // The result of the calculation will be stored here.

    std::vector<scalar> f(
        rgrid.size() );  // where the complete "potential" will be stored
    std::vector<scalar> wf( rgrid.size() );  // where the wf will be stored

    it<scalar> history;  // the last/"max" iteration info:
    it<scalar> current;


    current.iteration = 0;
    // current.excited = false;
    // current.upper_converged = false;
    // initialize to NaN's.... makes sure that things are done correctly
    for ( auto& i : wf ) i = std::numeric_limits<scalar>::quiet_NaN();
    for ( auto& i : f ) i  = std::numeric_limits<scalar>::quiet_NaN();

    // the last value of the potentialis the highest the energy can get.
    // This breaks down for continuum states
    current.energy_upper = pot( rgrid[wf.size() - 1], state );

    // the energy_upper isn't converged yet:
    current.upper_converged = false;
    current.energy_lower    = 0;
    for ( size_t i = 0; i < f.size(); i++ ) {
        current.energy_lower =
            std::min( current.energy_lower,
                      std::pow( scalar( state.l ) + .5, 2 ) /
                              ( 2. * std::pow( rgrid[i], 2 ) ) +
                          pot( rgrid[i], state ) );
    }
    // but the lowest the energy can get is the ground state or the last energy
    // level.
    // std::cerr << " energy_lower guess: " << current.energy_lower <<
    // std::endl;
    // std::cerr << " gs_energy guess: " << state.e.real() << std::endl;
    current.energy_lower =
        std::max( scalar( state.e.real() ), current.energy_lower );

    if ( current.energy_upper < current.energy_lower )
        current.energy_upper = current.energy_lower + 1;
    // The energy guess is an average of the lowest and the highest, but biased
    // towards the highest:
    current.energy = ( 5 * current.energy_upper + current.energy_lower ) / 11;

    // Then we screw things up:
    if ( current.energy_upper - current.energy_lower < 1. )
        current.energy_upper = current.energy_lower + 1.;
    // else if (current.energy_upper - current.energy_lower > 1.)
    // current.energy_upper = 2;

    current.de = 1e-3;  // we don't want to converge to quickly
    std::cerr << state << std::endl;
    std::cerr << "energy_upper: " << current.energy_upper
              << " energy_lower: " << current.energy_lower
              << " energy: " << current.energy << std::endl;

    bool converged = false;
    history        = current;
    history.energy = state.e.real();
    history.nodes  = -1;

    while ( !converged && current.iteration < 10000 ) {
        std::cerr << current << std::endl;
        current.iteration++;
        f[0] = 1 +
               dx * dx / 12 *
                   ( -std::pow( ( static_cast<scalar>( state.l ) + .5 ), 2 ) -
                     2 * std::pow( rgrid[0], 2 ) *
                         ( pot( rgrid[0], state ) - current.energy ) );

        std::vector<std::array<int, 2>> wells;
        std::array<int, 2>              w = {0, -1};
        bool in_well = ( f[0] - 1 < 0 ) ? false : true;
        std::cerr << "start in well? " << ( in_well ? "true" : "false" )
                  << std::endl;
        for ( size_t i = 1; i < f.size(); i++ ) {
            f[i] = dx * dx / 12 *
                   ( -std::pow( ( static_cast<scalar>( state.l ) + .5 ), 2 ) -
                     2 * std::pow( rgrid[i], 2 ) *
                         ( pot( rgrid[i], state ) - current.energy ) );

            if ( ( f[i] < 0 && f[i - 1] - 1 >= 0 ) )  // going out of well:
            {
                if ( in_well )  // I should be...
                {
                    w[1] = i;
                    wells.push_back( w );
                    in_well = false;
                } else
                    std::cerr << "I thought I was in well... but I'm not"
                              << std::endl;
            } else if ( ( f[i] > 0 &&
                          f[i - 1] - 1 <= 0 ) )  // going into a well
            {
                if ( in_well ) {
                    std::cerr << "I thought I was out of a well, but I'm not"
                              << std::endl;
                }
                if ( !in_well ) {
                    w       = std::array<int, 2>{{int( i ), -1}};
                    in_well = true;
                }
            }
            f[i] += 1;
        }
        if ( in_well == true ) {
            w[1] = rgrid.size() - 1;
            wells.push_back( w );
        }

        std::cerr << "wells: " << std::endl;
        for ( auto a : wells ) {
            std::cerr << a[0] << ", " << a[1] << std::endl;
            std::cerr << rgrid[a[0]] << ", " << rgrid[a[1]] << std::endl;
        }
        if ( wells.size() == 0 ) {
            w[1] = rgrid.size();
            wells.push_back( w );
        }

        if ( wells.size() >= 2 && state.j == 2 * state.l - 1 ) {
            messiness = wells[0][1];
        } else if ( wells.size() == 1 && state.j == 2 * state.l - 1 ) {
            if ( rgrid[wells[0][1]] < 1. )
                messiness = wells[0][1];
            else
                messiness = wells[0][0] / 2;
        } else
            messiness = 10;

        w = *std::max_element(
            wells.begin(), wells.end(),
            [=]( std::array<int, 2> a, std::array<int, 2> b ) {
                return ( ( rgrid[a[1]] - rgrid[a[0]] ) <
                         ( rgrid[b[1]] - rgrid[b[0]] ) );
            } );

        current.turnover = w[1];

        if ( w[1] != -1 && w[1] < rgrid.size() )
            std::cerr << " turnover: " << current.turnover
                      << " r[turnover]: " << rgrid[current.turnover]
                      << std::endl;
        else
            std::cerr << " turnover: " << current.turnover << " wells.size() "
                      << wells.size() << " wells: " << wells.front()[0] << ", "
                      << wells.front()[1] << std::endl;

        // The wavefunction should NOT have a classical turning point before
        // this:
        if ( current.turnover <= messiness ||
             current.turnover > ( rgrid.size() - 2 ) ) {
            if ( current.turnover != -1 && current.turnover < 2 ) {
                history              = current;
                current.energy_upper = current.energy;
                current.energy =
                    ( current.energy_upper + current.energy_lower ) / 2;
                current.nodes = -1;
                continue;
            } else if ( current.turnover > rgrid.size() - 2 ) {
                std::cerr << "excited!" << std::endl;
                // current.excited = true;
                // This is now an excited state...
                // current.energy_upper += std::abs( current.energy -
                // current.energy_lower );
            } else if ( current.turnover == -1 ||
                        current.turnover <= messiness ) {
                std::cerr << "the messiness is getting in the way" << std::endl;
                history              = current;
                current.energy_lower = current.energy;
                current.energy =
                    ( current.energy_upper + current.energy_lower ) / 2;
                current.nodes = -1;
                continue;
            }
        }
        // Set the initial values:

        // if the state is not a "problematic" state:
        if ( state.j != 2 * state.l - 1 || state.j == 0 ) {
            wf[0] = std::pow( rgrid[0], state.l + 1 ) *
                    ( 1. - 2. * rgrid[0] / ( 2. * scalar( state.l ) + 2. ) ) /
                    std::sqrt( rgrid[0] );
            wf[1] = std::pow( rgrid[1], state.l + 1 ) *
                    ( 1. - 2. * rgrid[1] / ( 2. * scalar( state.l ) + 2. ) ) /
                    std::sqrt( rgrid[1] );
            try {
                current.nodes = numerov(
                    f.begin(), ( f.begin() + current.turnover + 2 > f.end() ?
                                     f.end() :
                                     f.begin() + current.turnover + 2 ),
                    wf.begin(), ( wf.begin() + current.turnover + 2 > wf.end() ?
                                      wf.end() :
                                      wf.begin() + current.turnover + 2 ) );
            } catch ( std::out_of_range e ) {
            }

        } else  // if the state is:
        {
            wf[messiness] =
                std::pow( rgrid[messiness], state.l + 1 ) *
                ( 1. -
                  2. * rgrid[messiness] / ( 2. * scalar( state.l ) + 2. ) ) /
                std::sqrt( rgrid[messiness] );
            wf[messiness + 1] = std::pow( rgrid[messiness + 1], state.l + 1 ) *
                                ( 1. -
                                  2. * rgrid[messiness + 1] /
                                      ( 2. * scalar( state.l ) + 2. ) ) /
                                std::sqrt( rgrid[messiness + 1] );
            // zero everything before our new start
            for ( size_t i = 0; i < messiness; i++ ) wf[i] = 0;
            current.nodes =
                numerov( f.begin() + messiness,
                         ( f.begin() + current.turnover + 2 > f.end() ?
                               f.end() :
                               f.begin() + current.turnover + 2 ),
                         wf.begin() + messiness,
                         ( wf.begin() + current.turnover + 2 > wf.end() ?
                               wf.end() :
                               wf.begin() + current.turnover + 2 ) );
        }

        std::cerr << "nodes: " << current.nodes << std::endl;
        // now we check the nodes and iterate:
        if ( current.nodes - ( state.n - state.l - 1 ) >= 1 &&
             !current.upper_converged ) {
            std::cerr << "we are too high in energy: " << std::endl;
            if ( std::abs( current.de ) < 1e-18 && current.de == current.de ) {
                std::cerr << "subtract by de" << std::endl;
                current.energy -= std::abs( current.de );
                continue;
            }
            // we are too high in energy:
            if ( history.nodes == current.nodes - 1 && history.nodes != -1 ) {
                std::cerr << "set energy_upper to energy and average again, "
                             "cause we flipped"
                          << std::endl;
                auto t               = current;
                current.energy_upper = current.energy;
                current.de =
                    ( current.energy_upper + 3 * history.energy ) / 4 -
                    current.energy_upper;  // be biased towards the history end
                current.energy += current.de;
                history = t;
                continue;
            }
            std::cerr << "average and try again" << std::endl;

            history              = current;
            current.energy_upper = current.energy;
            // average the energy_upper and energy_lower and try again:
            current.de = ( current.energy_upper + current.energy_lower ) / 2 -
                         current.energy_upper;
            current.energy += current.de;
            continue;
        } else if ( current.nodes - ( state.n - state.l - 1 ) >= 1 &&
                    current.upper_converged ) {
            current.energy_lower = current.energy;
            current.energy =
                ( current.energy_upper + current.energy_lower ) / 2;
            continue;
        } else if ( current.nodes - ( state.n - state.l - 1 ) == 0 &&
                    !current.upper_converged ) {
            std::cerr << "we have the right energy" << std::endl;
            // current.upper_converged is false, so we actually want to go back
            // up:
            history = current;
            // difference between this and the last:
            // current.de = (current.energy + current.energy_upper)/2 -
            // current.energy;

            // if the "turnover point" is at the end of the grid, look at the
            // current.de.
            // we are now using that as our "end point":
            // if we are very close to the last upper limit that we had, we have
            // converged
            // otherwise, add the de and bisect:
            if ( std::abs( current.de ) < 1e-18 &&
                 current.turnover > wf.size() - 2 )  // if we have converged for
                                                     // an excited state:
            {
                converged = true;
                // figure out when the sign changes:
                if ( wf[4] < 0 )  // flip the wf:
                    for ( auto& a : wf ) a = -a;

                std::cout << current.turnover << "\t";
                continue;
            } else if ( std::abs( current.de ) < 1e-10 &&
                        current.turnover < wf.size() - 2 )  // we have converged
            // for a bound state,
            // time to do more
            // iterations:
            {
                std::cerr << "we have converged!";
                std::cerr << " de: " << current.de << std::endl;
                current.upper_converged = true;
                current.energy_upper    = current.energy;
                // we need a starting de, this shouldn't be changed:
                current.de = std::abs(
                    ( current.energy_upper - current.energy_lower ) / 2 );
                // This will continue below:
            } else if ( current.turnover > wf.size() - 2 ) {
                current.energy_lower = current.energy;
                current.de = ( current.energy_upper + current.energy ) / 2 -
                             current.energy;
                current.energy += current.de;
                continue;
            } else {
                current.de = ( current.energy_upper + current.energy ) / 2 -
                             current.energy;
                current.energy += current.de;
                // if (current.excited)
                // current.energy_lower = current.energy;
                continue;
            }
        }
        // if we are below (should only happen for lowest n for an l state):
        else if ( current.nodes - ( state.n - state.l - 1 ) <= -1 ) {
            std::cerr << "we are too low in energy: " << std::endl;
            if ( current.turnover > wf.size() - 2 )  // we are probably looking
                                                     // at an excited state, and
                                                     // we can't get to higher
                                                     // energy because of the
                                                     // current upper bound.
            {
                // bump the current energy upper bound by the de:
                history              = current;
                current.energy_lower = current.energy;
                current.energy =
                    ( current.energy_lower + current.energy_upper ) / 2;
                continue;
            }
            history              = current;
            current.energy_lower = current.energy;
            current.energy =
                ( current.energy_upper + current.energy_lower ) / 2;
            // if upper converged is true, then we have "moved too low", need to
            // readjust the de as well:
            if ( current.upper_converged )
                current.de = std::abs(
                    ( current.energy_upper - current.energy_lower ) / 2 );

            continue;
        }
        // if the state is converged... (shouldn't get here otherwise
        if ( !current.upper_converged )
            std::cerr << "how did I get here?" << std::endl;

        // start the end:
        wf[wf.size() - 1] = dx;
        wf[wf.size() - 2] = ( 12 - f[f.size() - 1] * 10 ) * wf[wf.size() - 1] /
                            f[wf.size() - 2];

        // save the points around the matching point:
        std::copy( wf.begin() + current.turnover - 1,
                   wf.begin() + current.turnover + 2, fore.begin() );
        // and figure out the rescaling factor:
        scalar rescale = wf[current.turnover];

        // run backwards from the end:
        numerov( f.rbegin(), f.rend() - current.turnover + 1, wf.rbegin(),
                 wf.rend() - current.turnover + 1 );
        std::copy( wf.begin() + current.turnover - 1,
                   wf.begin() + current.turnover + 2, back.begin() );

        // and make sure that the rescaling point is the same for both, and the
        // point before it is from
        // the forward iteration:
        wf[current.turnover - 1] = fore[0];

        // rescale:
        rescale = wf[current.turnover] / rescale;
        for ( size_t i = current.turnover; i < wf.size(); i++ )
            wf[i] /= rescale;
        for ( auto& a : back ) a /= rescale;

        // normalize:
        scalar norm = 0;
        for ( size_t i = 0; i < wf.size(); i++ )
            norm += wf[i] * wf[i] * rgrid[i] * rgrid[i] * dx;
        norm = std::sqrt( norm );
        for ( auto& a : wf ) a /= norm;
        for ( auto& a : fore ) a /= norm;
        for ( auto& a : back ) a /= norm;

        // find the difference in the derivatives:
        scalar forward_deriv1  = ( fore[0] - fore[1] ) / ( 2 * dx );
        scalar forward_deriv2  = ( fore[1] - fore[2] ) / ( 2 * dx );
        scalar backward_deriv1 = ( back[0] - back[1] ) / ( 2 * dx );
        scalar backward_deriv2 = ( back[1] - back[2] ) / ( 2 * dx );
        current.f              = -( forward_deriv1 + forward_deriv2 ) / 2 +
                    ( backward_deriv1 + backward_deriv2 ) / 2;

        // if the last iteration DIDN'T have the upper_converged, then we need
        // to set the "history"
        // to be the current iteration number:
        if ( history.upper_converged == false ) {
            history = current;
            // new energy will be halfway between the upper and lower:
            current.energy =
                ( current.energy_upper + current.energy_lower ) / 2;
            // current.de = std::abs( (current.energy - current.energy_upper)/2
            // );
            continue;
        }

        // otherwise, we do our convergence check:
        if ( std::abs( current.f ) < err || current.iteration > 600 ) {
            converged = true;
            std::cout << current.turnover << "\t";
            continue;
        }
        // then iterate if we haven't converged:
        else if ( history.f * current.f > 0 ) {
            history              = current;
            current.energy_upper = current.energy;
            current.energy =
                ( current.energy_lower + current.energy_upper ) / 2;
            continue;
        } else {
            current.energy_lower = current.energy;
            current.energy =
                ( current.energy_lower + current.energy_upper ) / 2;
            continue;
        }
    }
    if ( !converged ) {
        std::cerr << "didn't converge, returning anyways" << std::endl;
        std::cout << current.turnover << "*" << current.de << " ";
    }

    for ( size_t i = 0; i < wf.size(); i++ ) {
        wf[i] *= std::sqrt( rgrid[i] );
    }

    result.wf     = wf;
    result.energy = current.energy;
    return result;
};

template <typename scalar, typename write_type>
std::vector<BasisID>
n_loop( std::promise<std::complex<double>>&& future_guess, BasisID tmp,
        const std::vector<scalar>& rgrid,
        const BasisParameters<scalar, write_type>&      params,
        const std::function<scalar( scalar, BasisID )>& pot, scalar dx )
{
    std::vector<BasisID> energies( 0 );
    basis<scalar>        res;

    for ( tmp.n = tmp.l + 1; tmp.n <= params.nmax(); tmp.n++ ) {
        // if (params.bound_only() && temp.e.real() > 0.)
        // break;
        if ( params.fs() )
            for ( tmp.j = ( ( tmp.l > 0 ) ? 2 * tmp.l - 1 : 1 );
                  tmp.j <= ( ( tmp.l > 0 ) ? 2 * tmp.l + 1 : 1 ); tmp.j += 2 ) {
                std::cout << tmp << ", ";
                res   = find_basis<scalar>( tmp, dx, rgrid, pot, 1e-13 );
                tmp.e = res.energy;  // the energy min for the next will be the
                                     // correct energy for the last.
                // if (params.bound_only() && temp.e.real() > 0.)
                // break;
                energies.push_back( tmp );
                std::cout << ",\t" << res.energy << std::endl;
                // we need to convert the wf to PetscReal, or PetscScalar...
                math::normalize( res.wf, rgrid );
                common::export_vector_binary(
                    params.basis_function_filename( tmp ),
                    common::vector_type_change<scalar, write_type>( res.wf ) );
                std::cerr << tmp << std::endl;
                std::cerr << "========================================="
                          << std::endl;
                std::cerr << "========================================="
                          << std::endl;
                std::cerr << "========================================="
                          << std::endl;
                std::cerr << "========================================="
                          << std::endl;
                std::cerr.flush();
                if ( tmp.n == tmp.l + 2 &&
                     tmp.j == ( ( tmp.l > 0 ) ? 2 * tmp.l - 1 : 1 ) ) {
                    std::cout << "[" << std::this_thread::get_id()
                              << "] sending future " << tmp.l << ", " << tmp.j
                              << std::endl;
                    try {
                        future_guess.set_value( tmp.e );
                    } catch ( const std::future_error& e ) {
                        std::cout << "future error caught: " << e.code()
                                  << std::endl
                                  << e.what() << std::endl;
                    }
                    std::cout << "[" << std::this_thread::get_id()
                              << "] sent future" << std::endl;
                }
            }
        else {
            tmp.j = 0;
            std::cout << tmp << ", ";
            res   = find_basis<scalar>( tmp, dx, rgrid, pot, 1e-13 );
            tmp.e = res.energy;  // the energy min for the next will be the
                                 // correct energy for the last.
            energies.push_back( tmp );
            std::cout << ",\t" << res.energy << std::endl;
            // we need to convert the wf to PetscReal, or PetscScalar...
            math::normalize( res.wf, rgrid );
            common::export_vector_binary(
                params.basis_function_filename( tmp ),
                common::vector_type_change<scalar, write_type>( res.wf ) );
            std::cerr << tmp << std::endl;
            std::cerr << "========================================="
                      << std::endl;
            std::cerr << "========================================="
                      << std::endl;
            std::cerr << "========================================="
                      << std::endl;
            std::cerr << "========================================="
                      << std::endl;
            if ( tmp.n == tmp.l + 2 ) {
                std::cout << "[" << std::this_thread::get_id()
                          << "] sending future: " << tmp.l << std::endl;
                try {
                    future_guess.set_value( tmp.e );
                } catch ( const std::future_error& e ) {
                    std::cout << "future error caught: " << e.code()
                              << std::endl
                              << e.what() << std::endl;
                }
                std::cout << "[" << std::this_thread::get_id()
                          << "] sent future" << std::endl;
            }
        }
    }
    return energies;
}

template <typename scalar, typename write_type>
void find_basis_set( std::function<scalar( scalar, BasisID )> pot,
                     BasisParameters<scalar, write_type>&     params,
                     sae<scalar> atom )
{
    int rank, num;
    MPI_Comm_rank( params.comm(), &rank );
    MPI_Comm_size( params.comm(), &num );

    // Make grid:
    scalar              xmin = std::log( params.rmin() );
    scalar              xmax = std::log( params.rmax() );
    std::vector<scalar> xgrid( params.points() );
    for ( size_t i = 0; i < xgrid.size(); i++ )
        xgrid[i]   = xmin + i * ( xmax - xmin ) / ( params.points() - 1 );
    scalar dx      = xgrid[1] - xgrid[0];

    // get rgrid vector pointer from parameters and screw with it.
    std::vector<scalar>* rgrid = params.grid();
    for ( size_t i = 0; i < rgrid->size(); i++ ) {
        rgrid->at( i ) = std::exp( xgrid[i] );
    }

    std::vector<BasisID>* energies = params.basis_prototype();

    params.save_parameters();
    basis<scalar> res;
    BasisID       tmp;
    tmp.n = 1;
    tmp.l = 0;
    tmp.m = 0;
    tmp.j = 1;
    tmp.e = atom.gs_energy - 1;

    // concurrency and promises:
    size_t num_threads = params.procs();
    // std::promise<BasisID> p;
    // std::future<BasisID> f = p.get_future();
    // auto f2 = std::async( n_loop<scalar, write_type>, std::ref(p), tmp,
    // std::cref(*rgrid), std::cref(params), std::cref(pot), dx );

    energies->resize( 0 );
    //*energies = std::move( f2.get() );


    std::cout << "Number of threads: " << num_threads << std::endl;
    if ( rank == 0 ) std::cout << "n\tj\tl\te" << std::endl;
    std::cout << std::scientific;

    // std::cout << f.get() << std::endl;

    // tmp.e = f.get().e;
    for ( size_t l = 0; l <= params.lmax();
          l += num_threads )  // we will do some loop unrolling to make things
                              // simpler:
    {
        std::vector<std::future<std::vector<BasisID>>> futures_que(
            num_threads );
        // for each run through the "n's", start with an energy min of the gs,
        for ( size_t i = 0;
              i < ( l + num_threads > params.lmax() ? params.lmax() - l + 1 :
                                                      num_threads );
              i++ ) {
            std::promise<std::complex<double>> p_loop;
            std::future<std::complex<double>>  f_loop = p_loop.get_future();
            tmp.l                                     = l + i;
            futures_que[i] =
                std::async( std::launch::async, n_loop<scalar, write_type>,
                            std::move( p_loop ), tmp, std::cref( *rgrid ),
                            std::cref( params ), std::cref( pot ), dx );
            if ( tmp.l < params.nmax() - 1 ) {
                std::cout << "[0] waiting for future: " << tmp.l << std::endl;
                try {
                    f_loop.wait();
                    tmp.e = f_loop.get();
                } catch ( const std::future_error& e ) {
                    std::cout << "future error caught: " << e.code()
                              << std::endl
                              << e.what() << std::endl;
                }
                std::cout << "[0] got future" << std::endl;
            } else
                break;
        }
        std::cout << " after " << std::endl;
        for ( auto& a : futures_que ) {
            std::vector<BasisID> b;
            if ( a.valid() ) try {
                    a.wait();
                    b = a.get();
                } catch ( const std::future_error& e ) {
                    std::cout << "future error caught: " << e.code()
                              << std::endl
                              << e.what() << std::endl;
                }
            else
                continue;
            std::cout << " got future " << std::endl;
            energies->insert( energies->end(), b.begin(), b.end() );
        }
    }
    std::sort( energies->begin(), energies->end() );
    for ( auto& a : *energies ) std::cout << a << std::endl;
};
}

// We want to paralize this:
// we need:
// * a threadsafe "energies" container
// * a threadpool (so we aren't using more cores than we have)
// * a different way of structuring the for loop:  we want to start a number of
// "l" threads.  but only after we have already started the "l-1" thread
//
// int num_threads = std::thread::hardware_concurrency();
// std::vector< std::future< std::vector<BasisID> > > futures_que{ num_threads
// };
// std::promise< BasisID > p;
// tmp.e = atom.gs_energy - 1;
// f = p.get_future();
// futures_que[i] = n_loop( f, tmp, rgrid, params, pot, dx );


// for (size_t l = 1; l <= params.lmax(); l += num_threads) //we will do some
// loop unrolling to make things simpler:
//{
////for each run through the "n's", start with an energy min of the gs,
// std::future<BasisID> f;
// for (size_t i = 0; i < num_threads; i++ )
//{
// tmp.l = l + i;
// if (l + i == 0)
//{
// tmp.e = atom.gs_energy - 1;
// f = p.get_future();
// futures_que[i] = n_loop( f, tmp, rgrid, params, pot, dx );
//}
// else //or the l == l-1 state energy we have already found
//{

//}
//}


//};
