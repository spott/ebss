#pragma once

#include<paramters.hpp>

template <typename Scalar>
class Basis

{
protected:
    std::vector<Scalar> wf;

}


template <typename Scalar>
class Numerov: public Basis<Scalar>
{

    find_basis( Scalar (*pot)(Scalar r), BasisParameters params, int l)
    {
    }



private:
    int iterate(
}
