#pragma once

#include<vector>
#include<iterator>
#include<cmath>
#include<iostream>
//#include<common/parameters/Parameters.hpp>
#include<common/parameters/BasisParameters.hpp>
#include<common/common.hpp>
#include<limits>

namespace numerov
{

    enum type_node { A, B, C, D, E };

    template <typename iterator>
    int numerov(iterator fstart, iterator fend, iterator wfstart, iterator wfend)
    {
        int nodes = 0;

        if (wfend - wfstart != fend - fstart || wfend-wfstart == 0)
            std::cerr << "wf and f don't have the same size or are zero." << std::endl;

        bool inf = false;
        for(int i = 2; i < wfend-wfstart; i++)
        {
            wfstart[i] = ((12. - fstart[i-1] * 10.) * wfstart[i-1] - fstart[i-2] * wfstart[i-2]) / fstart[i];

            //std::cerr << i << " " << wfstart[i] << " ";
            if (*(wfstart+i) < 0 && *(wfstart+i-1) >= 0)
                nodes++;
            if (*(wfstart+i) > 0 && *(wfstart+i-1) <= 0)
                nodes++;
            //check for inf & nan:
            if (std::abs(*(wfstart+i)) == std::numeric_limits<typename iterator::value_type>::infinity() || (wfstart[i]) != (wfstart[i]))
            {
                std::cerr << "infinity detected at " << i << " last two values: " << *(wfstart+i-1) << ", " << wfstart[i-2] << std::endl;
                
                inf = true;
                break;
            }
        }

        return nodes;
    };

    template<typename scalar>
    struct basis {
        std::vector<scalar> wf;
        scalar energy;
    };

    template <typename scalar>
    basis<scalar> find_continuum(const int n, const int l, const scalar dx, const std::vector<scalar> &rgrid, scalar (*pot)(scalar r), bool &converged)
    {
        scalar energy_upper = 10;
        scalar energy_lower = -1;
        scalar energy = (energy_upper - energy_lower) / 2;
        int nodes;

        scalar de = 1;
        scalar err = 10e-17;

        std::vector<scalar> f(rgrid.size());
        std::vector<scalar> wf(rgrid.size());
        std::vector<scalar> potential(rgrid.size());

        //initialize to NaNs:
        for(size_t i = 0; i < wf.size(); i++)
            wf[i] = std::numeric_limits<scalar>::quiet_NaN();

        //initialize pot
        for(size_t i = 0; i < potential.size(); i++)
            potential[i] = pot(rgrid[i]);

        int iterations = -1;
        converged = false;

        type_node w;

        scalar norm;
        std::cerr.precision(22);
        while ( !converged )
        {
            iterations++;
            std::cerr << std::scientific <<  "excited: energy_upper: " << energy_upper << " energy_lower: " << energy_lower << " energy " << energy << std::endl ;

            //initialize f for the energy we are using
            for (size_t i = 0; i < f.size(); i++)
                f[i] = 1 + dx * dx / 12 * ( - std::pow((static_cast<scalar>(l) + .5), 2)
                        - 2 * std::pow(rgrid[i],2) * (pot(rgrid[i]) - energy));

            //set the initial values
            wf[0] = std::pow(rgrid[0],l+1) * (1. - 2. * rgrid[0] / (2. * (double)l + 2. )) / std::sqrt(rgrid[0]);
            wf[1] = std::pow(rgrid[1],l+1) * (1. - 2. * rgrid[1] / (2. * (double)l + 2. )) / std::sqrt(rgrid[1]);

            nodes = numerov(f.begin(),f.end(),wf.begin(),wf.end());

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
                //std::cerr << w << std::endl;
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
                //std::cerr << w << std::endl;
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
        int nan = -1;
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

    template <typename scalar>
    basis<scalar> find_bound(const int n, const int l, const scalar dx, const std::vector<scalar> &rgrid, scalar (*pot)(scalar r), bool &converged)
    {
        scalar de = 10e-10;
        scalar de2 = 10e-3;
        scalar err = 10e-20;
        scalar rescale, deriv;
        scalar norm;

        std::vector<scalar> f(rgrid.size());
        std::vector<scalar> wf(rgrid.size());
        std::vector<scalar> potential(rgrid.size());

        //initialize to NaNs:
        for(size_t i = 0; i < wf.size(); i++)
            wf[i] = std::numeric_limits<scalar>::quiet_NaN();

        //initialize pot
        //for(size_t i = 0; i < potential.size(); i++)
            //potential[i] = pot(rgrid[i]);

        scalar energy_upper = pot(rgrid[rgrid.size()-1]);
        scalar energy_lower=0;
        scalar energy = 0;
        for (size_t i = 1; i < rgrid.size(); i++)
        {
            energy_lower = std::min(energy_lower, (std::pow((static_cast<scalar>(l)+.5),2)/(2. * std::pow(rgrid[i],2)) + pot(rgrid[i])));
        }
        energy = (energy_upper + energy_lower) / 2;


        int nodes;
        int turnover = -1;
        int iterations = -1;
        converged = false;
        int deriverror;

        std::cerr.precision(22);
        while ( !converged )
        {
            iterations++;
            std::cerr << std::scientific <<  "energy_upper: " << energy_upper << " energy_lower: " << energy_lower << " energy " << energy << " de: " << de << " de2: " << de2 << std::endl ;

            //initialize f for the energy we are using
            f[0] = 1 + dx * dx / 12 * ( - std::pow((static_cast<scalar>(l) + .5), 2)
                    - 2 * std::pow(rgrid[0],2) * (pot(rgrid[0]) - energy));
            for (size_t i = 1; i < f.size(); i++)
            {
                f[i] = dx * dx / 12 * ( - std::pow((static_cast<scalar>(l) + .5), 2)
                        - 2 * std::pow(rgrid[i],2) * (pot(rgrid[i]) - energy));
                if (f[i] > 0 && f[i-1] - 1 <= 0)
                    turnover = i;
                if (f[i] < 0 && f[i-1] - 1 >= 0)
                    turnover = i;
                f[i] += 1;
            }

            if (turnover < 2 || turnover > static_cast<int>(rgrid.size()) - 9)
            {
                energy -= de2 + de;
                de2 /= 2;
                std::cerr << "oops (" << iterations << ")" << std::endl;
                if (iterations >= 1000)
                {
                    break;
                }
                continue;
            }

            //set the initial values
            wf[0] = std::pow(rgrid[0],l+1) * (1. - 2. * rgrid[0] / (2. * (double)l + 2. )) / std::sqrt(rgrid[0]);
            wf[1] = std::pow(rgrid[1],l+1) * (1. - 2. * rgrid[1] / (2. * (double)l + 2. )) / std::sqrt(rgrid[1]);

            //std::cerr << f[3] << " " << wf[1] << std::endl;

            wf[wf.size()-1] = dx;
            wf[wf.size()-2] = (12 - f[f.size()-1] * 10) * wf[wf.size()-1] / f[wf.size()-2];
            
            //std::cerr << wf[wf.size() - 1] << " " << wf[wf.size() - 2] << std::endl;
            nodes = numerov(f.begin(), f.begin() + turnover + 1, wf.begin(), wf.begin() + turnover + 1);

            //std::cerr << "wf: " << wf[0] << " " << wf[3] << " " << wf[turnover] <<std::endl;
            rescale = wf[turnover];

            numerov(f.rbegin(), f.rend() - turnover , wf.rbegin(), wf.rend() - turnover );
            
            //std::cerr << "wf: " << wf[turnover-2] << " " << wf[turnover-1] << " " << wf[turnover] <<std::endl;

            rescale = wf[turnover]/rescale;

            for (size_t i = turnover; i < wf.size(); i++)
                wf[i] /= rescale;

            //normalize!
            norm = 0.0;
            for (size_t i = 0; i < wf.size(); i++)
                norm += wf[i] * wf[i] * rgrid[i] * rgrid[i] * dx;
            norm = std::sqrt(norm);
            for (size_t i = 0; i < wf.size(); i++)
                wf[i] /= norm;

            scalar cusp = (wf[turnover-1] * f[turnover-1] + 
                    f[turnover+1] * wf[turnover+1] + 
                    10 * f[turnover] * wf[turnover]) / 12;
            scalar dfcusp = f[turnover] * (wf[turnover] / cusp - 1);
            de = dfcusp * 12 * dx * std::pow(cusp,2);


            if (nodes != n - l - 1 )
            {
                //std::cerr << w << std::endl;
                de2 /= 2;
                if (nodes > n - l - 1)
                    energy_upper = energy;
                if (nodes < n - l - 1)
                    energy_lower = energy;
                energy = ( energy_upper + energy_lower )/2;
                if (iterations >= 10000)
                {
                    break;
                }
                continue;
            }

            deriv = (wf[wf.size()-1] - wf[wf.size()-2]) / dx;
            //stderr.writeln("deriv: ",deriv);
            if (abs(deriv) > 10e-5)
            {
                deriverror++;
                if (deriverror >= 20)
                    break;
            }

            //std::cerr << std::scientific <<  "energy_upper: " << energy_upper << " energy_lower: " << energy_lower << " energy " << energy << " de: " << de << " de2: " << de2 << std::endl ;
            //std::cerr << "nodes: " << nodes << " norm: " << norm << " rescale " << rescale << " cusp: " << cusp << " dfcusp " << dfcusp << " deriv " << deriv << std::endl;
            if (de < 0 && de2 > 0)
                de2 = .5 * std::abs(de2) * de / std::abs(de);
            if (de > 0 && de2 < 0)
                de2 = .5 * std::abs(de2) * de / std::abs(de);

            if (de > 0.)
                energy_lower = energy;
            if (de < 0.)
                energy_upper = energy;

            if (std::abs(de) > std::abs(de2) )
                energy += de;
            if (std::abs(de2) > std::abs(de) )
                energy += de2;

            //check for convergence
            if (std::abs(de) < std::abs(err))
                converged = true;
            if (iterations > 5000)
                break;

            //make sure we don't go out of energy bounds
            energy = std::min(energy, energy_upper);
            energy = std::max(energy, energy_lower);
        }

        //multiply by sqrt(rgrid[]), required because of log grid:
        int nan = -1;
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

    int this_l_excited_n = 1000;
    int this_l = 1000;

    template <typename scalar>
    basis<scalar> find_basis(const int n, const int l, const scalar dx, const std::vector<scalar> &rgrid, scalar (*pot)(scalar r))
    {
        bool converged = false;
        basis<scalar> result;
        if (l != this_l)
        {
            this_l = l;
            this_l_excited_n = 1000;
        }
        if (n < this_l_excited_n)
        {
            result = find_bound(n, l, dx, rgrid, pot, converged);
            if (converged == false)
            {
                this_l_excited_n = n;
                //std::cout << "excited: " << std::endl;
                result = find_continuum(n, l, dx, rgrid, pot, converged);
            }
        }
        else
        {
            //std::cout << "excited: " << std::endl;
            result = find_continuum(n, l, dx, rgrid, pot, converged);
        }

        return result;
    }

    template<typename scalar, typename write_type>
    void find_basis_set( scalar (*pot)(scalar), BasisParameters<scalar, write_type> *params)
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

        params->save_parameters();
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
                std::vector<write_type> wf2 = common::vector_type_change<scalar, write_type>(res.wf);
                common::export_vector_binary(params->basis_function_filename(n,l), &wf2); 
            }

        //Need to combine the energy vectors... and ideally sort them...

    };


}
