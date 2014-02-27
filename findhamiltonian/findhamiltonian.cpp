
//ebss:
#include<common/parameters/HamiltonianParameters.hpp>
#include<common/common.hpp>
#include<common/math.hpp>

//petsc:
#include<petsc.h>

//stl:
#include<vector>
#include<sstream>
#include<string>
#include<iostream>
#include<unordered_map>
#include<memory>
#include<functional>

//gsl:
#include<gsl/gsl_sf_coulomb.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_errno.h>


template <typename ReturnType, typename Arg>
std::function<ReturnType (Arg)> memoize(std::function<ReturnType (Arg)> func)
{
    std::unordered_map<Arg, ReturnType> cache;
    return ([=](Arg t) mutable  {
            if (cache.size() > 400)
                cache.clear();
            if (cache.find(t) == cache.end())
                cache.emplace(t,func(t));
            return cache[t];
    });
}

template <typename ReturnType, typename Arg>
std::function<ReturnType (Arg, bool)> half_memoize(std::function<ReturnType (Arg)> func, size_t size)
{
    std::unordered_map<Arg, ReturnType> cache;
    return ([=](Arg t, bool remember) mutable  {
            if (cache.find(t) == cache.end())
            {
                //if (remember == false)
                    //return func(t);
                if (cache.size() > size)
                {
                    if (remember)
                    {
                        cache.clear();
                        cache.emplace(t,func(t));
                    }
                    else
                        return func(t);
                }
                if (cache.size() > size / 2 && !remember)
                    return func(t);
                cache.emplace(t,func(t));
            }
            return cache[t];
            });
}

struct dipole_params { BasisID a,b; };

double dipole_matrix_element( double r, void* params )
{
    dipole_params* p = static_cast<dipole_params*>(params);

//std::cout <<  p->a << " -- " << p->b << " -- " << r << std::endl;

    gsl_sf_result left;
    gsl_sf_result right;
    auto error_handler = gsl_set_error_handler_off();
    //auto left = gsl_sf_hydrogenicR (p->a.n, p->a.l, 1, r);
    //auto right = gsl_sf_hydrogenicR (p->b.n, p->b.l, 1, r);
    auto err = gsl_sf_hydrogenicR_e(p->a.n, p->a.l, 1, r, &left);
    if (err == GSL_EUNDRFLW)
        left.val = 0;
    else if (err != 0) {
        std::cout << err << " " << gsl_strerror(err) << std::endl;
        throw std::exception();
    }

    err = 0;
    err = gsl_sf_hydrogenicR_e(p->b.n, p->b.l, 1, r, &right);
    if (err == GSL_EUNDRFLW)
        right.val = 0;
    else if (err != 0) {
        std::cout << err << " " << gsl_strerror(err) << std::endl;
        throw std::exception();
    }

    gsl_set_error_handler (error_handler);

    return left.val * r * r * r * right.val;
}

int main(int argc, const char ** argv)
{
    int ac = argc;
    char** av = new char*[argc];
    for (size_t i = 0; i < argc; i++)
    {
        av[i] = new char[strlen(argv[i])];
        std::copy(argv[i], argv[i] + strlen(argv[i]), av[i]);
    }
    PetscInitialize(&ac, &av, PETSC_NULL, PETSC_NULL);

    HamiltonianParameters<PetscReal> params(argc, argv, MPI_COMM_WORLD);
    if (params.rank() == 0)
        std::cout << params.print();

    //decltype( params.basis_parameters()->grid() ) grid;
    std::vector<double> grid;
    auto prototype = params.prototype();

    if (! params.analytic() )  {
        grid = params.basis_parameters()->grid();
    }

    // we need to do a rough calculation of how much memory we are going to need to hold, in order
    // to figure out how large to make our memoization routine:
    
    //total size of dipole_matrix:
    auto memory = prototype.size()*(params.fs() ? 2 : 1) * 2 * params.nmax() * sizeof(PetscScalar);
    if (params.rank() == 0) std::cout << "total size of dipole matrix: " << memory;

    //size per proc (fudge factor of 2 so we don't get too big):
    memory /= 2 * params.size();
    if (params.rank() == 0) std::cout << ", per processor: " << memory;

    //we have ~2gb per proc:
    memory = params.mem_per_proc() - memory;
    if (params.rank() == 0) std::cout << ", left over per processor for memoization: " << memory;
    
    //each
    size_t per_basis_state = 0;
    if (! params.analytic() )
    {
        per_basis_state = grid.size() * sizeof(double) + 100;
        if (params.rank() == 0) std::cout << std::endl << "each basis state takes up: " << per_basis_state ;
    }
    else
        per_basis_state = 1;

    //left over number of basis states:
    memory /= per_basis_state;

    if (params.rank() == 0) std::cout << ", leaving enough for: " << memory << " basis states in memoization" <<std::endl;

    std::function<bool (int, int)> dipole_selection_rules;
    if (params.fs())
    {
        dipole_selection_rules = [prototype](int i, int j) {
            // since j is multiplied by 2, need to have a 2 between them.
            return ( (std::abs(prototype[i].l - prototype[j].l) == 1 && 
                        (std::abs(prototype[i].j - prototype[j].j) == 2 || 
                        prototype[i].j == prototype[j].j )) || 
                    i == j);
        };
    } 
    else
    {
        dipole_selection_rules = [prototype](int i, int j) {
            return (std::abs(prototype[i].l - prototype[j].l) == 1 || i == j );
        };
    }

    std::function< std::shared_ptr<std::vector<double> > ( BasisID ) > import_wf = 
        [&params](BasisID a) -> std::shared_ptr<std::vector<double> > {
            std::shared_ptr< std::vector<double> > b(
                    new std::vector<double>(std::move(
                            common::import_vector_binary<double>(
                                params.basis_parameters()->basis_function_filename(a)
                                )
                            )
                        )); 
            return b;
        };
    auto import_wf2 = half_memoize(import_wf, memory);


    std::function<PetscScalar (int, int)> findvalue;
    if (params.fs())
    {
        findvalue = [&prototype,&grid,&import_wf2](int i, int j)->PetscScalar{
            if (i == j)
                return 0.0;
            auto a = import_wf2(prototype[i], true);
            auto b = import_wf2(prototype[j], false);
            int s = 0;
            for (size_t n = 0; n < (*a).size(); ++n)
            {
                if (math::signum( (*a)[n] ) != 0 && math::signum( (*b)[n] ) != 0)
                {
                    s =  math::signum( (*a)[n] * (*b)[n] );
                    break;
                }
            }

            //if (s == 0)  std:: cout << "shit... one or the other wavefunction is always zero" << std::endl;
            PetscScalar radial = math::integrateTrapezoidRule(
                    *a,
                    *b,
                    grid );
            PetscScalar angular = math::CGCoefficient<PetscScalar>(prototype[i],prototype[j]);
            //we are only considering spin up electrons.  so m_j always == +1/2 (since m_l is always 0)
            angular *= std::sqrt (
                    (prototype[i].l + (prototype[i].j - 2 * prototype[i].l) * .5 * .5 + .5) * 
                    (prototype[j].l + (prototype[j].j - 2 * prototype[j].l) * .5 * .5 + .5) / 
                    ((2. * prototype[i].l + 1) * (2. * prototype[j].l + 1)));
            if (prototype[i].n == 2 && prototype[j].n == 2)
                std::cout << "n = 2 transition: " << radial * angular << std::endl;
            return radial * angular;
        };
    }
    else if ( !params.analytic() )
    {
        findvalue = [&prototype,&grid,&import_wf2](int i, int j)->PetscScalar{
            if (i == j)
                return 0.0;
            auto a = import_wf2(prototype[i], true);
            auto b = import_wf2(prototype[j], false);

            int s = 0;
            for (size_t n = 0; n < (*a).size(); ++n)
            {
                if (std::abs( (*a)[n] ) > 1e-6 && std::abs( (*b)[n] ) >  1e-6)
                {
                    s =  math::signum( (*a)[n] * (*b)[n]);
                    //std::cout << s << " from: " << (*a)[n] << " and " << (*b)[n] << std::endl;
                    break;
                }
            }

            //if (s == 0)  std:: cout << "shit... one or the other wavefunction is always zero" << std::endl;
            PetscScalar radial = math::integrateTrapezoidRule(
                    *a,
                    *b,
                    grid );
            PetscScalar angular = math::CGCoefficient<PetscScalar>(prototype[i],prototype[j]);
            if (prototype[i].n == 2 && prototype[j].n == 2)
                std::cout << "n = 2 transition: " << radial * angular << std::endl;
            return radial * angular;
        };
    }
    else
    {
        findvalue = [&prototype](int i, int j)->PetscScalar{
            if (i == j)
                return 0.0;
            dipole_params p{prototype[i], prototype[j]};
            double radial = 0;
            double error;
            gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
            gsl_function f;
            f.function = &dipole_matrix_element;
            f.params = static_cast<void*>( &p );

            auto error_handler = gsl_set_error_handler_off();
            std::cout << p.a << " --- " << p.b << "\t";
            auto err = gsl_integration_qagiu (&f, 0, 1e-8, 1e-9, 10000, w, &radial, &error);
            std::cout << radial << "\t" << error << "\t|| " << gsl_strerror(err) << std::endl;
            gsl_set_error_handler (error_handler);

            PetscScalar angular = math::CGCoefficient<PetscScalar>(prototype[i],prototype[j]);
            if (prototype[i].n == 2 && prototype[j].n == 2)
                std::cout << "n = 2 transition: " << radial * angular + 3 << std::endl;
            gsl_integration_workspace_free(w);
            return radial * angular;
        };


    }


   Mat H = common::populate_matrix< std::complex<double> >(params, 
                                    dipole_selection_rules, 
                                    findvalue, 
                                    prototype.size(),
                                    prototype.size(), true);


   params.write_dipole_matrix(H);
   params.write_energy_eigenvalues();

   params.save_parameters();
   //delete params;
   PetscFinalize();
   return 0;
}
