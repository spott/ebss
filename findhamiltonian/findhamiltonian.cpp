
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

    auto grid = params.basis_parameters()->grid();
    auto prototype = params.prototype();

    // we need to do a rough calculation of how much memory we are going to need to hold, in order
    // to figure out how large to make our memoization routine:
    
    //total size of dipole_matrix:
    auto memory = prototype.size()*(params.fs() ? 2 : 1) * 2 * params.nmax() * sizeof(PetscScalar);
    if (params.rank() == 0) std::cout << "total size of dipole matrix: " << memory;

    //size per proc (fudge factor of 2 so we don't get too big):
    memory /= 2 * params.size();
    if (params.rank() == 0) std::cout << ", per processor: " << memory;

    //we have ~2gb per proc:
    memory = 193274000 - memory;
    if (params.rank() == 0) std::cout << ", left over per processor for memoization: " << memory;
    
    //each
    auto per_basis_state = grid->size() * sizeof(double) + 100;
    if (params.rank() == 0) std::cout << std::endl << "each basis state takes up: " << per_basis_state ;

    //left over number of basis states:
    memory /= per_basis_state;

    if (params.rank() == 0) std::cout << ", leaving enough for: " << memory << " basis states in memoization" <<std::endl;


    std::function<bool (int, int)> dipole_selection_rules;
    if (params.fs())
    {
        dipole_selection_rules = [prototype](int i, int j) {
            // since j is multiplied by 2, need to have a 2 between them.
            return ( std::abs(prototype[i].l - prototype[j].l) == 1 && (
                        std::abs(prototype[i].j - prototype[j].j) == 2 || 
                        prototype[i].j == prototype[j].j ) || 
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
            PetscScalar radial = std::abs(math::integrateTrapezoidRule(
                    *a,
                    *b,
                    *grid ));
            PetscScalar angular = math::CGCoefficient<PetscScalar>(prototype[i],prototype[j]);
            //we are only considering spin up electrons.  so m_j always == +1/2 (since m_l is always 0)
            angular *= std::sqrt ( 
                    (prototype[i].l + (prototype[i].j - 2 * prototype[i].l) * .5 * .5 + .5) * 
                    (prototype[j].l + (prototype[j].j - 2 * prototype[j].l) * .5 * .5 + .5) / 
                    ((2. * prototype[i].l + 1) * (2. * prototype[j].l + 1)));
            return radial * angular;
        };
    }
    else
    {
        findvalue = [&prototype,&grid,&import_wf2](int i, int j)->PetscScalar{
            if (i == j)
                return 0.0;
            auto a = import_wf2(prototype[i], true);
            auto b = import_wf2(prototype[j], false);
            PetscScalar radial = math::integrateTrapezoidRule(
                    *a,
                    *b,
                    *grid );
            PetscScalar angular = math::CGCoefficient<PetscScalar>(prototype[i],prototype[j]);
            return radial * angular;
        };
    }


   Mat H = common::populate_matrix< std::complex<double> >(params, 
                                    dipole_selection_rules, 
                                    findvalue, 
                                    prototype.size(),
                                    prototype.size());


   params.write_dipole_matrix(H);
   params.write_energy_eigenvalues();

   params.save_parameters();
   //delete params;
   PetscFinalize();
   return 0;
}
