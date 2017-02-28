#pragma once

#include <cmath>
#include <common/common.hpp>
#include <common/math.hpp>
#include <common/parameters/BasisParameters.hpp>
#include <functional>
#include <unordered_map>

template <typename T>
struct sae_param {
    int p;
    T   c, beta;
};

template <typename T>
struct sae {
    std::vector<sae_param<T>> params;
    double                    Z;
    double                    N;
    T                         gs_energy;
};


template <typename T>
T parameterized_pot( const T r, const sae<T>& atom, const BasisID& state )
{
    T a = 0;
    for ( auto p : atom.params )
        a += p.c * std::pow( r, p.p ) * std::exp( -p.beta * r );
    a *= ( atom.N - 1 );
    a += 1 - atom.N + atom.Z;
    a *= -1 / r;
    return a;
}

template <typename T>
T parameterized_finestructure_pot( const T r, const sae<T>& atom,
                                   const BasisID& state )
{
    T a = parameterized_pot( r, atom, state );

    T b = 0;
    for ( auto p : atom.params )
        b += p.c * std::pow( r, p.p - 2 ) * std::exp( -p.beta * r ) *
             ( 1. - p.p + r * p.beta );
    b *= T( atom.N - 1 );
    b += 1. / ( r * r );
    b *= T( ( T( state.j ) / 2. ) * ( T( state.j ) / 2. + 1. ) -
            T( state.l ) * ( T( state.l ) + 1. ) - 3. / 4. );
    b *= 1. / ( 4. * math::C * math::C * r );
    return a + b;
}

template <typename T>
T helium_tf( const T r, const BasisID& state )
{
    T a = 1.0 + 1.231 * std::exp( -0.662 * r ) -
          1.325 * r * std::exp( -1.236 * r ) - 0.231 * std::exp( -0.480 * r );
    a *= -1.0;
    a /= r;
    return a;
}

template <typename T>
T argon_tf( const T r, const BasisID& state )
{
    T a = 1.0 + 16.039 * std::exp( -2.007 * r ) -
          25.543 * r * std::exp( -4.525 * r ) + 0.961 * std::exp( -0.443 * r );
    a *= -1.0;
    a /= r;
    return a;
}

template <typename T>
T krypton_michelle( const T r, const BasisID& state )
{
    // Michelle SAE potential:
    // V(r) = -c0./r - (Zc.*exp(-r0.*r))./r - a1.*exp(-b1.*r) - a2.*exp(-b2.*r)
    // - a3.*exp(-b3.*r) - a4.*exp(-b4.*r)
    // 1	35	1.78886177	-984.162936	18.613213	951.023767	18.958577	-61.9173561
    // 3.29525736	3.49098308	199.999977
    T a = -1.0/r - 35. * std::exp( -1.78886177 * r )/r +
          984.162936 * std::exp( -18.613213 * r ) -
          951.023767 * std::exp( -18.958577 * r ) +
          61.9173561 * std::exp( -3.29525736 * r ) -
          3.49098308 * std::exp( -199.999977 * r );
    return a;
}


template <typename T>
std::function<T( const T, BasisID )> memoized_pot( const sae<T>& atom )
{
    // BasisID state;
    // auto cache = std::make_shared<std::unordered_map< T, T> >();
    auto func = std::bind( parameterized_pot<T>, std::placeholders::_1,
                           std::cref( atom ), std::placeholders::_2 );
    return func;
    // return ( [state, [>cache,<] func](const T r, BasisID s) mutable {
    // return func(r, state);
    ////if (s != state)
    ////{
    ////cache->clear();
    ////state = s;
    ////(*cache)[r] = func(r, s);
    ////return (*cache)[r];
    ////}
    ////else
    ////{
    ////if (cache->find(r) == cache->end())
    ////{
    ////(*cache)[r] = func(r, state);
    ////}
    ////return (*cache)[r];
    ////}
    //});
}

template <typename T>
std::function<T( const T, BasisID )>
memoized_finestructure_pot( const sae<T>& atom )
{
    // BasisID state;
    // auto cache = std::make_shared<std::unordered_map< T, T> >();
    auto func =
        std::bind( parameterized_finestructure_pot<T>, std::placeholders::_1,
                   std::cref( atom ), std::placeholders::_2 );
    return func;
    // return ( [state, [>cache,<] func](const T r, BasisID s) mutable {
    // return func(r, state);
    // if (s != state)
    //{
    // cache->clear();
    // state = s;
    //(*cache)[r] = func(r, s);
    // return (*cache)[r];
    //}
    // else
    //{
    // if (cache->find(r) == cache->end())
    //{
    //(*cache)[r] = func(r, state);
    //}
    // return (*cache)[r];
    //}
    //});
}
