#pragma once

#include<vector>
#include<array>
#include<iterator>
#include<cmath>
#include<iostream>
#include<common/common.hpp>
#include<limits>
#include<common/math.hpp>
#include<common/output.hpp>

template <typename T>
struct sae_param {
    int p;
    T c, beta;
};

template <typename T>
struct sae {
    std::vector < sae_param<T> > params;
    int Z;
    int N;
    T gs_energy;
};

namespace numerov
{


    template <typename iterator>
    int numerov(iterator fstart, iterator fend, iterator wfstart, iterator wfend)
    {
        unsigned int nodes = 0;

        if (wfend - wfstart != fend - fstart || wfend-wfstart == 0)
            std::cerr << "wf and f don't have the same size or are zero." << std::endl;

        bool inf = false;
        for(int i = 2; i < wfend-wfstart; i++)
        {
            if (!inf)
                wfstart[i] = ((12. - fstart[i-1] * 10.) * wfstart[i-1] - fstart[i-2] * wfstart[i-2]) / fstart[i];
            else
                wfstart[i] = wfstart[i-1];

            //std::cerr << i << " " << wfstart[i] << " ";
            if (*(wfstart+i) < 0 && *(wfstart+i-1) >= 0)
                nodes++;
            if (*(wfstart+i) > 0 && *(wfstart+i-1) <= 0)
                nodes++;
            //check for inf & nan:
            if (std::abs(*(wfstart+i)) == std::numeric_limits<typename iterator::value_type>::infinity() || (wfstart[i]) != (wfstart[i]))
            {
                wfstart[i] = wfstart[i-1];
                std::cerr << "infinity detected at " << i << " last two values: " << *(wfstart+i-1) << ", " << wfstart[i-2] << " nodes: " << nodes << std::endl;
                //throw (std::exception());
                inf = true;
            }
        }

        return nodes;
    };

    template <typename scalar>
        struct basis {
            std::vector<scalar> wf;
            scalar energy;
        };

    template <typename scalar>
        basis<scalar> 
        find_basis( const BasisID state, 
                    const scalar dx, const std::vector<scalar> &rgrid, 
                    const std::function< scalar (scalar, BasisID) > pot, 
                    const scalar err,
                    const scalar e_guess,
                    bool &converged )
        {
            struct it {
                int iteration;
                scalar energy;
                int nodes;
                int turnover;
                //std::array<int,2> well;
                // test values:
                std::array< scalar, 3 > fore;
                std::array< scalar, 3 > back;
                scalar deriv;
                scalar de;
                scalar f;
            };

            basis<scalar> result;

            std::vector<scalar> f(rgrid.size());
            std::vector<scalar> wf(rgrid.size());
            std::vector<scalar> potential(rgrid.size());

            //initialize to NaNs:
            for(size_t i = 0; i < wf.size(); i++)
                wf[i] = std::numeric_limits<scalar>::quiet_NaN();

            scalar energy_upper = pot(rgrid[rgrid.size()-1], state);
            scalar energy_lower=0;
            scalar energy = 0;
            for (size_t i = 1; i < rgrid.size(); i++)
            {
                energy_lower = std::min(energy_lower, (std::pow((static_cast<scalar>(state.l)+.5),2)/(2. * std::pow(rgrid[i],2)) + pot(rgrid[i], state)));
            }

            //write out potential:
            //std::vector<scalar> potential(rgrid->size());
            for (size_t i = 0; i < rgrid.size(); i++)
            {
                potential[i] = pot(rgrid[i], state);
            }
            std::stringstream ss;
            ss << "potential_l_" << state.l << "_j_" << state.j << ".dat";
            common::export_vector_ascii(ss.str(),potential);
            //if the energy lower bound is much much lower than the gs_energy, cut it off...
            if (energy_lower < state.e.real() - 200)
                energy_lower = state.e.real() - 200;

            energy = state.e.real();
                //(energy_lower + energy_upper) / 2;

            it now = {-1, energy, -1, -1, {0,0,0}, {0,0,0}, 0, 1., 10};
            it last = now;
            scalar rescale;
            scalar norm;
            converged = false;
            scalar fmax;

            while ( !converged )
            {
                now.iteration++;

                std::cerr << output::red << std::scientific <<  "energy_upper: " << energy_upper << " energy_lower: " << energy_lower << " energy " << now.energy << " deriv: " << now.deriv  << " iterations: " << now.iteration << " node: " << now.nodes << output::reset << std::endl;

                //initialize f:
                f[0] = 1 + dx * dx / 12 * ( - std::pow((static_cast<scalar>(state.l) + .5), 2)
                        - 2 * std::pow(rgrid[0],2) * (pot(rgrid[0], state) - now.energy));

                //find the classic "wells" (all of the paired turnover points)
                std::vector< std::array<int,2> > wells;
                std::array< int, 2 > w = {0,-1};
                bool in_well = true;
                for (size_t i = 1; i < f.size(); i++)
                {
                    f[i] = dx * dx / 12 * ( - std::pow((static_cast<scalar>(state.l) + .5), 2)
                            - 2 * std::pow(rgrid[i],2) * (pot(rgrid[i], state) - now.energy));

                    if (( (f[i] < 0 && f[i-1] - 1 >= 0) || (f[i] > 0 && f[i-1] -1 <= 0) ))
                    {
                        if (in_well)
                        {
                            w[1] = i;
                            wells.push_back(w);
                            in_well = false;
                        }
                        if (!in_well)
                        {
                            w = std::array<int,2>{int(i), -1};
                            in_well = true;
                        }
                    }
                    f[i] += 1;
                }
                
                std::cerr << "wells: " << std::endl;
                for (auto a: wells)
                    std::cerr << a[0] << ", " << a[1] << std::endl;


                if (wells.size() == 0)
                    wells.push_back(w);
                //figure out which well is biggest, and use the end of that as our "turnaround point".
                w = *std::max_element(wells.begin(),
                                 wells.end(), 
                                 [](std::array< int, 2 > a, std::array< int, 2 > b) 
                                            { return (a[1]-a[2] < b[1]-b[0]); } 
                                );
                now.turnover = w[1];
                if (w[1] != -1)
                    std::cerr << " turnover: " << now.turnover << " r[turnover]: " << rgrid[now.turnover] << std::endl;
                //TODO check the last iteration, and go to continuum states if they are both (or more?) 
                //failing
                if (now.turnover < 3000 || now.turnover > (rgrid.size() - 2) )
                {
                    if (now.turnover != -1 && now.turnover < 2)
                    {
                        energy_upper = now.energy;
                        last = now;
                        now.energy = (energy_upper + energy_lower)/2;
                        now.nodes = -1;
                        continue;
                    }
                    else if (now.turnover == -1 || now.turnover < 3000)
                    {
                        energy_lower = now.energy;
                        last = now;
                        now.energy = (energy_upper + energy_lower)/2;
                        now.nodes = -1;
                        continue;
                    }
                }

                //set the initial values
                if (state.j != 2 * state.l - 1)
                {
                    wf[0] = std::pow(rgrid[0],state.l+1) * (1. - 2. * rgrid[0] / (2. * scalar(state.l) + 2. )) 
                        / std::sqrt(rgrid[0]);
                    wf[1] = std::pow(rgrid[1],state.l+1) * (1. - 2. * rgrid[1] / (2. * scalar(state.l) + 2. )) 
                        / std::sqrt(rgrid[1]);
                    now.nodes = numerov(f.begin(), f.begin() + now.turnover + 2, wf.begin(), wf.begin() + now.turnover + 2);
                }
                else
                {
                    wf[3000] = std::pow(rgrid[3000],state.l+1) * (1. - 2. * rgrid[3000] / (2. * scalar(state.l) + 2. )) 
                        / std::sqrt(rgrid[3000]);
                    wf[3001] = std::pow(rgrid[3001],state.l+1) * (1. - 2. * rgrid[3001] / (2. * scalar(state.l) + 2. )) 
                        / std::sqrt(rgrid[3001]);
                    for (size_t i = 0; i < 3000; i++)
                        wf[i] = 0;
                    now.nodes = numerov(f.begin() + 3000, f.begin() + now.turnover + 2, wf.begin() + 3000, wf.begin() + now.turnover + 2);
                }
                wf[wf.size()-1] = dx;
                wf[wf.size()-2] = (12 - f[f.size()-1] * 10) * wf[wf.size()-1] / f[wf.size()-2];

                std::cerr << "nodes: " << now.nodes << std::endl;
                if (now.nodes > state.n - state.l - 1 )
                {
                    energy_upper = now.energy;
                    if (last.nodes <= state.n - state.l - 1)
                    {
                        std::cerr << " energy " <<  now.energy << " last.energy: " << last.energy << std::endl;
                        it t = now;
                        now.energy = (now.energy + last.energy)/2;
                        last = t;
                        if (std::abs(now.energy - last.energy) > 0.00)
                            now.de /= 2;
                        else
                            now.de *= .5;
                        now.nodes = -1;
                        std::cerr << " energy " <<  now.energy << " last.energy: " << last.energy << std::endl;
                    }
                    else
                    {
                        last = now;
                        now.nodes = -1;
                        now.energy = (energy_upper + energy_lower)/2;
                    }
                    std::cerr << "continuuing.... " << std::endl;
                    continue;
                }
                else if (now.nodes < state.n - state.l - 1 )
                {
                    energy_lower = now.energy;
                    if (last.nodes >= state.n - state.l - 1)
                    {
                        it t = now;
                        now.energy = (now.energy + last.energy)/2;
                        last = t;
                        now.de /= 2;
                        now.nodes = -1;
                    }
                    else
                    {
                        last = now;
                        now.nodes = -1;
                        now.energy = (energy_lower + energy_upper)/2 ;
                    }
                    //now.energy = (energy_upper + energy_lower)/2;
                    continue;
                }
                std::copy(wf.begin() + now.turnover - 1, wf.begin() + now.turnover + 2, now.fore.begin());
                //std::cerr << " fore: " << now.fore[0] << ", " << now.fore[1] << ", " << now.fore[2];
                //std::cerr << " wf: " << wf[now.turnover-1] << ", " << wf[now.turnover] << ", " << wf[now.turnover+1];
                rescale = wf[now.turnover];
                numerov(f.rbegin(), f.rend() - now.turnover + 1 , wf.rbegin(), wf.rend() - now.turnover + 1 );
                std::copy(wf.begin() + now.turnover - 1, wf.begin() + now.turnover + 2, now.back.begin());
                wf[now.turnover-1] = now.fore[0];

                std::cerr << " fore: " << now.fore[0] << ", " << now.fore[1] << ", " << now.fore[2];
                std::cerr << " back: " << now.back[0] << ", " << now.back[1] << ", " << now.back[2];
                std::cerr << std::endl;
                std::cerr << " wf: " << wf[now.turnover - 5] << ", "<< wf[now.turnover - 4] << ", "<< wf[now.turnover - 3] << ", "<< wf[now.turnover - 2] << ", "<< wf[now.turnover - 1] << ", *"<< wf[now.turnover] << ", "<< wf[now.turnover + 1] << ", "<< wf[now.turnover + 2] << ", "<< wf[now.turnover + 3] << ", "<< wf[now.turnover + 4] << ", "<< wf[now.turnover + 5] << ", " << std::endl;
                rescale = wf[now.turnover]/rescale;

                for (size_t i = now.turnover; i < wf.size(); i++)
                    wf[i] /= rescale;

                for (auto &a: now.back)
                    a /= rescale;

                norm = 0;
                for (size_t i = 0; i < wf.size(); i++)
                    norm += wf[i] * wf[i] * rgrid[i] * rgrid[i] * dx;
                norm = std::sqrt(norm);
                std::cerr << " norm: " << norm << std::endl;
                for (auto &a: wf)
                    a /= norm;
                for (auto &a: now.fore)
                    a /= norm;
                for (auto &a: now.back)
                    a /= norm;


                std::cerr << " fore: " << now.fore[0] << ", " << now.fore[1] << ", " << now.fore[2];
                std::cerr << " back: " << now.back[0] << ", " << now.back[1] << ", " << now.back[2];
                std::cerr << std::endl;

                scalar derivf1 = (now.fore[0] - now.fore[1]) /(2*dx);
                scalar derivf2 = (now.fore[1] - now.fore[2]) /(2*dx);
                scalar derivb1 = (now.back[0] - now.back[1]) /(2*dx);
                scalar derivb2 = (now.back[1] - now.back[2]) /(2*dx);
                if ( math::signum(derivf1) != math::signum(derivb2) )
                {
                    std::cerr << output::blue << "deriv1: " << derivf1 << ", " << math::signum(derivf1) << " deriv2 " << derivb2 << ", " << math::signum(derivb2) << output::reset <<std::endl;
                    last = now;
                    now.energy = now.energy + now.de;
                    now.energy = std::min(now.energy, energy_upper);
                    now.energy = std::max(now.energy, energy_lower);
                    continue;
                }


                scalar cusp = (wf[now.turnover-1] * f[now.turnover-1] + 
                               f[now.turnover+1] * wf[now.turnover+1] + 
                          10 * f[now.turnover] * wf[now.turnover]) / 12;
                scalar dfcusp = f[now.turnover] * (wf[now.turnover] / cusp - 1);
                now.deriv = dfcusp * 12 * dx * std::pow(cusp,2);
                std::cerr << "deriv: " <<  now.deriv << std::endl;
                std::cerr << "f: " << derivf1 - derivb2 << std::endl;
                std::cerr << "f2: " << derivb1 - derivf2 << std::endl;
                std::cerr << "f3: " << (derivf1 + derivf2)/2 - (derivb1 + derivb2)/2 << std::endl;
                now.f = (derivf1 + derivf2)/2 - (derivb1 + derivb2)/2;


                //flip back and forth, but keep the de big:
                if (math::signum(last.deriv) == - math::signum(now.deriv))
                {
                    if (now.energy > last.energy)
                        energy_lower = last.energy;
                    else if (now.energy < last.energy)
                    {
                        energy_upper = last.energy; 
                        fmax = last.f;
                    }
                    std::cerr << "now.energy - last.energy " << now.energy - last.energy << std::endl;
                    if (std::abs(now.energy - last.energy) >= 0.000) 
                        now.de = (now.energy + last.energy)/2 - now.energy;
                }
                else if (fmax * now.f > 0)
                {
                    if (std::abs((derivf1 + derivf2)/2 - (derivb1 + derivb2)/2) < 1e-15)
                    {
                        std::cerr << output::green << "left because of convergence: de:" << now.de << " err: " << err << output::reset << std::endl;
                        converged = true;
                    }
                    last = now;
                    now.energy -= std::abs(now.de);
                    continue;
                }
                else
                {
                    if (std::abs((derivf1 + derivf2)/2 - (derivb1 + derivb2)/2) < 1e-14)
                    {
                        std::cerr << output::green << "left because of convergence: de:" << now.de << " err: " << err << output::reset << std::endl;
                        converged = true;
                    }
                    last = now;
                    now.de /= 2;
                    now.energy += std::abs(now.de);
                    continue;
                }
                //{
                    //std::cerr << " going the right direction, or wrong direction: " << last.deriv << ", " << now.deriv<< ", de:" << now.de  <<std::endl;
                    //now.energy += std::abs(now.de) * math::signum(now.deriv);
                    //now.de /= 2;
                //}


                if (std::abs((derivf1 + derivf2)/2 - (derivb1 + derivb2)/2) < 1e-15)
                {
                    std::cerr << output::green << "left because of convergence: de:" << now.de << " err: " << err << output::reset << std::endl;
                    converged = true;
                }


                last = now;
                //last.energy = std::min(last.energy, energy_upper);
                //last.energy = std::max(last.energy, energy_lower);

                now.energy += std::abs(now.de) * math::signum(now.deriv);

                now.energy = std::min(now.energy, energy_upper);
                now.energy = std::max(now.energy, energy_lower);
            }
            std::cout << now.turnover << "\t";
            result.wf = wf;
            result.energy = now.energy;
            return result;
        }

    template <typename scalar, typename write_type>
        void find_basis_set( std::function< scalar (scalar, BasisID) > pot, 
                             BasisParameters<scalar, write_type> *params, 
                             sae<scalar> atom)
        {
            int rank, num;
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

            std::vector<BasisID> *energies = params->basis_prototype();

            params->save_parameters();
            basis<scalar> res;
            BasisID tmp;

            energies->resize(0);
            //MPI stuff to split up workload:
            bool converged = false;

            if (rank==0) std::cout << "n\tl\tj\tturn\te" << std::endl;
            std::cout << std::scientific;
            for (int l = 0; l <= params->lmax(); l += num)
            {
                //int n = 4;
                for (int n = l+1; n <= params->nmax(); n++)
                {
                    for (int j = ((l>0)? 2*l+1 : 1 ); j >= ((l>0)? 2 * l - 1 : 1); j-=2)
                    {
                        //tmp.e = -.159;
                        if (1 == n)
                            tmp.e = atom.gs_energy;
                        else
                            tmp.e = atom.gs_energy;
                        tmp.n = n;
                        tmp.j = j;
                        tmp.m = 0;
                        tmp.l = l;
                        std::cout << tmp.n << "\t" << tmp.l << "\t" << tmp.j << "\t";
                        res = find_basis<scalar>( tmp, dx, *rgrid, pot, 1e-18, tmp.e.real(), converged);
                        tmp.e = res.energy;
                        energies->push_back(tmp);
                        std::cout << res.energy << std::endl;
                        //we need to convert the wf to PetscReal, or PetscScalar...
                        common::export_vector_binary(
                                params->basis_function_filename(tmp),
                                common::vector_type_change<scalar, write_type>(res.wf));
                        std::cerr << n << "\t" << l << "\t" << j << "\t" << res.energy << std::endl;
                        std::cerr << "=========================================" << std::endl;
                        std::cerr << "=========================================" << std::endl;
                        std::cerr << "=========================================" << std::endl;
                        std::cerr << "=========================================" << std::endl;
                    }
                }
            }
        };
}


