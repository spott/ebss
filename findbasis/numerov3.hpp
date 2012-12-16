#pragma once

#include<vector>
#include<forward_list>
#include<array>
#include<iterator>
#include<cmath>
#include<iostream>
#include<common/common.hpp>
#include<limits>
#include<common/math.hpp>

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
    struct it {
        int iteration;
        scalar energy_upper; //Eventually will be the n-l-1 == 1 energy
        bool upper_converged; //this tells us if the energy_upper is converged on the correct value;
        scalar energy_lower; //This won't change, except for special cases (but we will test for it anyways...
        scalar energy;      //the energy of this iteration
        int nodes;          //how many nodes for this iteration
        int turnover;       //where the turnover is...
        scalar f;       //the derivative matching value;
        scalar de;  //Will be used to bisect
    };


    template <typename scalar>
    std::ostream& operator<<(std::ostream &out, const it<scalar> &b)     //output
    {
        out << "iteration: " << b.iteration << " energy_upper: " << b.energy_upper << " upper_converged: " << (b.upper_converged ? "true" : "false")  << " energy_lower: " << b.energy_lower << "\n\t energy: " << b.energy << " f: " << b.f << " de: " << b.de;
        return out;
    };

    template <typename scalar>
        basis<scalar> 
        find_basis( const BasisID state, 
                    const scalar dx, const std::vector<scalar> &rgrid, 
                    const std::function< scalar (scalar, BasisID) > pot, 
                    const scalar err)
        {

            int messiness = 3000; //where the messiness for the shitty j's is adjusted away:
            std::array< scalar, 3 > fore;  //I don't think I need a memory of these...
            std::array< scalar, 3 > back;
            basis<scalar> result; //The result of the calculation will be stored here.

            std::vector<scalar> f(rgrid.size());  // where the complete "potential" will be stored
            std::vector<scalar> wf(rgrid.size()); // where the wf will be stored

            it<scalar> history;  //the last/"max" iteration info:
            it<scalar> current;


            current.iteration = 0;
            //current.upper_converged = false;
            //initialize to NaN's.... makes sure that things are done correctly
            for (auto &i: wf)
                i = std::numeric_limits<scalar>::quiet_NaN();
            for (auto &i: f)
                i = std::numeric_limits<scalar>::quiet_NaN();

            //the last value of the potentialis the highest the energy can get.  
            //This breaks down for continuum states
            current.energy_upper = pot(rgrid[ wf.size()-1 ],state);
            //the energy_upper isn't converged yet:
            current.upper_converged = false;
            current.energy_lower = 0;
            for (size_t i = 0; i < f.size(); i++)
            {
                current.energy_lower = std::min (
                        current.energy_lower,
                        std::pow( scalar(state.l) + .5, 2) / (2. * std::pow(rgrid[i],2) ) + pot(rgrid[i], state)
                        );
            }
            //but the lowest the energy can get is the ground state or the last energy level.
            current.energy_lower = std::max(scalar(state.e.real()), current.energy_lower);
            //The energy guess is an average of the lowest and the highest, but biased towards the highest:
            current.energy = (10 * current.energy_upper + current.energy_lower)/11;
            std::cerr << state << std::endl;
            std::cerr << "energy_upper: " << current.energy_upper << " energy_lower: " << current.energy_lower << " energy: " << current.energy << std::endl;

            bool converged = false;

            while (!converged)
            {
                std::cerr << current << std::endl;
                current.iteration++;
                f[0] = 1 + dx * dx / 12 * ( - std::pow((static_cast<scalar>(state.l) + .5), 2)
                        - 2 * std::pow(rgrid[0],2) * (pot(rgrid[0], state) - current.energy));

                std::vector< std::array<int,2> > wells;
                std::array< int, 2 > w = {0,-1};
                bool in_well = (f[0] - 1 < 0)? false : true;
                for (size_t i = 1; i < f.size(); i++)
                {
                    f[i] = dx * dx / 12 * ( - std::pow((static_cast<scalar>(state.l) + .5), 2)
                            - 2 * std::pow(rgrid[i],2) * (pot(rgrid[i], state) - current.energy));

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

                if (wells.size() >=2 && state.j == 2 * state.l - 1)
                    messiness = wells[0][1];
                else
                    messiness = 3000;

                w = *std::max_element(wells.begin(),
                                 wells.end(), 
                                 [](std::array< int, 2 > a, std::array< int, 2 > b) 
                                            { return (a[1]-a[2] < b[1]-b[0]); } 
                                );
                current.turnover = w[1];

                if (w[1] != -1)
                    std::cerr << " turnover: " << current.turnover << " r[turnover]: " << rgrid[current.turnover] << std::endl;

                //The wavefunction should NOT have a classical turning point before this:
                if (current.turnover < messiness || current.turnover > (rgrid.size() - 2) )
                {
                    if (current.turnover != -1 && current.turnover < 2)
                    {
                        history = current;
                        current.energy_upper = current.energy;
                        current.energy = (current.energy_upper + current.energy_lower)/2;
                        current.nodes = -1;
                        continue;
                    }
                    else if (current.turnover == -1 || current.turnover < messiness)
                    {
                        history = current;
                        current.energy_lower = current.energy;
                        current.energy = (current.energy_upper + current.energy_lower)/2;
                        current.nodes = -1;
                        continue;
                    }
                }

                //Set the initial values:

                //if the state is not a "problematic" state:
                if (state.j != 2 * state.l - 1)
                {
                    wf[0] = std::pow(rgrid[0],state.l+1) * (1. - 2. * rgrid[0] / (2. * scalar(state.l) + 2. )) 
                        / std::sqrt(rgrid[0]);
                    wf[1] = std::pow(rgrid[1],state.l+1) * (1. - 2. * rgrid[1] / (2. * scalar(state.l) + 2. )) 
                        / std::sqrt(rgrid[1]);
                    current.nodes = numerov(f.begin(), f.begin() + current.turnover + 2, wf.begin(), wf.begin() + current.turnover + 2);
                }
                else //if the state is:
                {
                    wf[messiness] = std::pow(rgrid[messiness],state.l+1) * (1. - 2. * rgrid[messiness] / (2. * scalar(state.l) + 2. )) 
                        / std::sqrt(rgrid[messiness]);
                    wf[messiness+1] = std::pow(rgrid[messiness+1],state.l+1) * (1. - 2. * rgrid[messiness+1] / (2. * scalar(state.l) + 2. )) 
                        / std::sqrt(rgrid[messiness+1]);
                    //zero everything before our new start
                    for (size_t i = 0; i < messiness; i++)
                        wf[i] = 0;
                    current.nodes = numerov(f.begin() + messiness, f.begin() + current.turnover + 2, wf.begin() + messiness, wf.begin() + current.turnover + 2);
                }

                std::cerr << "nodes: " << current.nodes << std::endl;
                //now we check the nodes and iterate:
                if (current.nodes - (state.n - state.l - 1) >= 1 && !current.upper_converged)
                {
                    //we are too high in energy:
                    history = current;
                    current.energy_upper = current.energy;
                    //average the energy_upper and energy_lower and try again:
                    current.de = (current.energy_upper + current.energy_lower)/2 - current.energy_upper;
                    current.energy += current.de;
                    continue;
                }
                else if (current.nodes - (state.n - state.l - 1) >= 1 && current.upper_converged)
                {
                    current.energy_lower = current.energy;
                    current.energy = (current.energy_upper + current.energy_lower)/2;
                    continue;
                }
                else if (current.nodes - (state.n - state.l - 1) == 0 && !current.upper_converged)
                {
                    //current.upper_converged is false, so we actually want to go back up:
                    history = current;
                    //difference between this and the last:
                    current.de = (current.energy + current.energy_upper)/2 - current.energy;
                    //if we are very close to the last upper limit that we had, we have converged
                    //otherwise, add the de and bisect:
                    if (current.de > 1e-5)
                    {
                        current.energy += current.de;
                        continue;
                    }
                    else
                    {
                        current.upper_converged = true;
                        current.energy_upper = current.energy;
                        //we need a starting de, this shouldn't be changed:
                        current.de =std::abs( (current.energy_upper - current.energy_lower) / 2);
                        //This will continue below:
                    }
                }
                // if we are below (should only happen for lowest n for an l state):
                else if (current.nodes - (state.n-state.l-1) <= -1)
                {
                    history = current;
                    current.energy_lower = current.energy;
                    current.energy = (current.energy_upper + current.energy_lower)/2;
                    //if upper converged is true, then we have "moved too low", need to readjust the de as well:
                    if (current.upper_converged)
                        current.de = std::abs( (current.energy_upper - current.energy_lower)/2);

                    continue;
                }
                //if the state is converged... (shouldn't get here otherwise
                if (!current.upper_converged)
                    std::cerr << "how did I get here?" << std::endl;

                //start the end:
                wf[wf.size()-1] = dx;
                wf[wf.size()-2] = (12 - f[f.size()-1] * 10) * wf[wf.size()-1] / f[wf.size()-2];

                //save the points around the matching point:
                std::copy(wf.begin() + current.turnover - 1, wf.begin() + current.turnover + 2, fore.begin());
                //and figure out the rescaling factor:
                scalar rescale = wf[current.turnover];

                //run backwards from the end:
                numerov(f.rbegin(), f.rend() - current.turnover + 1 , wf.rbegin(), wf.rend() - current.turnover + 1 );
                std::copy(wf.begin() + current.turnover - 1, wf.begin() + current.turnover + 2, back.begin());

                //and make sure that the rescaling point is the same for both, and the point before it is from 
                //the forward iteration:
                wf[current.turnover-1] = fore[0];

                //rescale:
                rescale = wf[current.turnover]/rescale;
                for (size_t i = current.turnover; i < wf.size(); i++)
                    wf[i] /= rescale;
                for (auto &a: back)
                    a /= rescale;

                //normalize:
                scalar norm = 0;
                for (size_t i = 0; i < wf.size(); i++)
                    norm += wf[i] * wf[i] * rgrid[i] * rgrid[i] * dx;
                norm = std::sqrt(norm);
                for (auto &a: wf)
                    a /= norm;
                for (auto &a: fore)
                    a /= norm;
                for (auto &a: back)
                    a /= norm;

                //find the difference in the derivatives:
                scalar forward_deriv1 = (fore[0] - fore[1]) /(2*dx);
                scalar forward_deriv2 = (fore[1] - fore[2]) /(2*dx);
                scalar backward_deriv1 = (back[0] - back[1]) /(2*dx);
                scalar backward_deriv2 = (back[1] - back[2]) /(2*dx);
                current.f = -(forward_deriv1 + forward_deriv2)/2 + (backward_deriv1 + backward_deriv2)/2;

                //if the last iteration DIDN'T have the upper_converged, then we need to set the "history"
                //to be the current iteration number:
                if (history.upper_converged == false)
                {
                    history = current;
                    //new energy will be halfway between the upper and lower:
                    current.energy = (current.energy_upper + current.energy_lower)/2;
                    //current.de = std::abs( (current.energy - current.energy_upper)/2 );
                    continue;
                }

                //otherwise, we do our convergence check:
                if (std::abs(current.f) < err || current.iteration > 500)
                {
                    converged = true;
                    std::cout << current.turnover << "\t";
                    continue;
                }
                //then iterate if we haven't converged:
                else if (history.f * current.f > 0)
                {
                    history = current;
                    current.energy_upper = current.energy;
                    current.energy = (current.energy_lower + current.energy_upper) / 2;
                    continue;
                }
                else
                {
                    current.energy_lower = current.energy;
                    current.energy = (current.energy_lower + current.energy_upper) / 2;
                    continue;
                }

            }
            result.wf = wf;
            result.energy = current.energy;
            return result;
        };

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

            if (rank==0) std::cout << "n\tl\tj\te" << std::endl;
            std::cout << std::scientific;
            for (int l = rank; l <= params->lmax(); l += num)
            {
                //for each run through the "n's", start with an energy min of the gs,
                if (l == 0)
                    tmp.e = atom.gs_energy - 1;
                else //or the l == 0 state energy we have already found
                {
                    tmp = *std::find_if(
                            energies->begin(), 
                            energies->end(), 
                            [l](const BasisID &a) {return (a.n == l+1 && a.l == l-1 && a.j == 2 * a.l + 1); }
                            );
                    tmp.e = tmp.e*1.1;
                }
                for (int n = l+1; n <= params->nmax(); n++)
                {
                    for (int j = ((l>0)? 2 * l - 1 : 1); j <= ((l>0)? 2*l+1 : 1 ); j+=2)
                    {
                        //tmp.e = -.159;
                        tmp.n = n;
                        tmp.j = j;
                        tmp.m = 0;
                        tmp.l = l;
                        std::cout << n << "\t" << l << "\t" << j << "\t";
                        res = find_basis<scalar>( tmp, dx, *rgrid, pot, 1e-13);
                        tmp.e = res.energy;     //the energy min for the next will be the correct energy for the last.
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
