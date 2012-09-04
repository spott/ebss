#pragma once

#include<vector>
#include<cmath>
#include<iostream>
#include<common/parameters.hpp>
#include<common/common.hpp>
#include<limits>

namespace numerov
{

    enum type_node { A, B, C, D, E };


    template <typename scalar>
    int numerov(std::vector<scalar> *f, std::vector<scalar> *wf)
    {
        int nodes = 0;

        if (wf->size() != f->size() || wf->size() == 0)
            std::cerr << "wf and f don't have the same size or are zero. wf: " 
                << (*wf).size() << ", f: " << (*f).size() << std::endl;

        for(size_t i = 2; i < wf->size(); i++)
        {
            (*wf)[i] = ((12. - (*f)[i-1] * 10.) * (*wf)[i-1] - (*f)[i-2] * (*wf)[i-2]) / (*f)[i];

            if ((*wf)[i] < 0 && (*wf)[i-1] >= 0)
                nodes++;
            if ((*wf)[i] > 0 && (*wf)[i-1] <= 0)
                nodes++;

        }

        return nodes;
    };
    template<typename scalar>
    struct basis {
        std::vector<scalar> wf;
        scalar energy;
    };

    template <typename scalar>
    basis<scalar> find_basis(const int n, const int l, const scalar dx, const std::vector<scalar> &rgrid, scalar (*pot)(scalar r))
    {
        scalar energy_upper = 10;
        scalar energy_lower = -1;
        scalar energy = (energy_upper - energy_lower) / 2;
        int nodes;

        scalar de = 1;
        scalar err = 10e-20;

        std::vector<scalar> f(rgrid.size());
        std::vector<scalar> wf(rgrid.size());
        std::vector<scalar> potential(rgrid.size());

        //initialize to NaNs:
        for(size_t i; i < wf.size(); i++)
            wf[i] = std::numeric_limits<scalar>::quiet_NaN();

        //initialize pot
        for(size_t i; i < potential.size(); i++)
            potential[i] = pot(rgrid[i]);

        int iterations = -1;
        bool converged = false;

        type_node w;

        scalar norm;
        std::cerr.precision(22);
        while ( !converged )
        {
            iterations++;
            //std::cerr << std::scientific <<  "energy_upper: " << energy_upper << " energy_lower: " << energy_lower << " energy " << energy << std::endl ;

            //initialize f for the energy we are using
            for (size_t i = 0; i < f.size(); i++)
                f[i] = 1 + dx * dx / 12 * ( - std::pow((static_cast<scalar>(l) + .5), 2) 
                        - 2 * std::pow(rgrid[i],2) * (pot(rgrid[i]) - energy));

            //set the initial values
            wf[0] = std::pow(rgrid[0],l+1) * (1. - 2. * rgrid[0] / (2. * (double)l + 2. )) / std::sqrt(rgrid[0]); 
            wf[1] = std::pow(rgrid[1],l+1) * (1. - 2. * rgrid[1] / (2. * (double)l + 2. )) / std::sqrt(rgrid[1]); 

            nodes = numerov<scalar>(&f, &wf);

            //normalize!
            norm = 0.0;
            for (size_t i = 0; i < wf.size(); i++)
                norm += wf[i] * wf[i] * rgrid[i] * rgrid[i] * dx;
            norm = std::sqrt(norm);
            for (size_t i = 0; i < wf.size(); i++)
                wf[i] /= norm;

            // Right: between A and B (just shy of the correct number of nodes +1)
            // more nodes: A, D
            // fewer nodes: B, C, E
            if (nodes == n - l )
            {
                std::cerr << w << std::endl;
                energy_upper = energy;
                if (w == C || w == E)
                    de /= -2;
                if (w == B)
                    de /= -2;

                energy = energy + de;
                w = A;
            }
            else if (nodes == n - l - 1)
            {
                std::cerr << w << std::endl;
                if (w == A)
                    de /= -2;
                if (w == D)
                    de *= -1;
                energy = energy + de;
                w = B;
            }
            else if (nodes == n - l - 2)
            {
                energy_lower = energy;
                if (w == A)
                    de /= -2;
                if (w == D)
                    de /= -2;
                energy += de;
                w = C;
            }
            else if (nodes > n - l)
            {
                if (w == A)
                    de *= -1;
                if (w == B || w == C)
                    de /= -2;
                if (w == D && de > 0)
                    de *= -1;
                if (w == E)
                    de /= -2;

                energy_upper = energy;
                energy += de;
                w = D;
            }
            else if (nodes < n - l - 2)
            {
                energy_lower = energy;

                if (w == E && de < 0)
                    de *= -1;
                else if (w == D || w == A)
                    de /= -2;
                energy += de;
                w = E;
            }
            //std::cerr << " de " << de << std::endl; 

            //check for convergence
            if (std::abs(de) < std::abs(err))
                converged = true;

            //make sure we don't go out of energy bounds
            energy = std::min(energy, energy_upper);
            energy = std::max(energy, energy_lower);
        }

        //multiply by sqrt(rgrid[]), required because of log grid:
        nan = -1;
        for (size_t i = 0; i < wf.size(); i++)
        {
            if (wf[i] != wf[i] || std::abs(wf[i]) == std::numeric_limits<scalar>::infinity())
                nan=i;
            wf[i] *= std::sqrt(rgrid[i]);
        }

        //Check for and report NaNs:
        if (nan != -1 )
            std::cerr << "the wavefunction is NaN... at: " << nan << " n: " << n << ", l: "<< l << ", e: " << energy << " sqrt "<<std::endl;

        //normalize!
        norm = 0.0;
        for (size_t i = 0; i < wf.size(); i++)
        {
            norm += wf[i] * wf[i] * rgrid[i] * dx;
            if (wf[i] != wf[i] || std::abs(wf[i]) == std::numeric_limits<scalar>::infinity())
                nan=i;
        }

        //check for nan's
        if (nan != -1 )
            std::cerr << "the wavefunction is NaN... at: " << nan << " n: " << n << ", l: "<< l << ", e: " << energy << " norm: " << norm << std::endl;

        norm = std::sqrt(norm);
        for (size_t i = 0; i < wf.size(); i++)
        {
            wf[i] /= norm;
            if (wf[i] != wf[i] || std::abs(wf[i]) == std::numeric_limits<scalar>::infinity())
                nan=i;
        }
        if (nan != -1 )
            std::cerr << "the wavefunction is NaN... at: " << nan << " n: " << n << ", l: "<< l << ", e: " << energy << std::endl;
        basis<scalar> out;
        out.wf = wf;
        out.energy = energy;

        return out;

    };


    template<typename scalar>
    void find_basis_set( scalar (*pot)(scalar), BasisParameters<scalar> *params)
    {
        int rank;
        int num;
        MPI_Comm_rank(params->comm(), &rank);
        MPI_Comm_size(params->comm(), &num);
        //Make grid:
        scalar xmin = std::log(params->rmin());
        scalar xmax = std::log(params->rmax());

        std::vector<scalar> xgrid(params->points());
        for (size_t i = 0; i < xgrid.size(); i++)
            xgrid[i] = xmin + i * (xmax-xmin)/params->points();

        scalar dx = xgrid[1] - xgrid[0];

        //get rgrid vector pointer from parameters and screw with it.
        std::vector<scalar> *rgrid = params->grid();
        for (size_t i = 0; i < rgrid->size(); i++)
        {
            rgrid->at(i) = std::exp(xgrid[i]);
        }


        //get energy vector pointer from parameters
        std::vector<BasisID> *energies = params->basis_prototype();

        basis<scalar> res;
        BasisID tmp;


        //MPI stuff to split up workload:

        if (rank==0) std::cout << "n\tl\te\t∆\t∆/exact" << std::endl;

        for (int l = rank; l <= params->lmax(); l += num)
            for (int n = l+1; n <= params->nmax(); n++)
            {
                res = find_basis(n, l, dx, *rgrid, pot );
                tmp.n = n;
                tmp.m = 0;
                tmp.l = l;
                tmp.e = res.energy;
                energies->push_back(tmp);
                std::cout << n << "\t" << l << "\t" << res.energy << "\t" << res.energy + 1./(2.*n*n) << "\t" << (res.energy + 1./(2.*n*n))/(1./(2.*n*n)) << std::endl;
                //we need to convert the wf to PetscReal, or PetscScalar...
                std::vector<PetscReal> wf2 = common::vector_type_change<scalar, PetscReal>(res.wf);
                common::export_vector_binary(params->basis_function_filename(n,l), &wf2); 
            }

        //Need to combine the energy vectors... and ideally sort them...

    };


}
