
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
            if (cache.size() > 200)
                cache.clear();
            if (cache.find(t) == cache.end())
                cache.emplace(t,func(t));
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

    HamiltonianParameters<PetscReal> *params = new HamiltonianParameters<PetscReal>(argc, argv, MPI_COMM_WORLD);
    if (params->rank() == 0)
        std::cout << params->print();

    std::vector<PetscReal> *grid = params->basis_parameters()->grid();
    auto prototype = params->prototype();

    std::function<bool (int, int)> dipole_selection_rules = [prototype,params](int i, int j) {
                if (params->with_fs())
                    // since j is multiplied by 2, need to have a 2 between them.
                    return ((std::abs(prototype[i].j - prototype[j].j) == 2) || prototype[i].j == prototype[j].j || i == j );
                else
                    return ((std::abs(prototype[i].l - prototype[j].l) == 1) || i == j );
    };

    std::function< std::shared_ptr<std::vector<double> > ( BasisID ) > import_wf = 
        [&params](BasisID a) -> std::shared_ptr<std::vector<double> > {
        std::shared_ptr< std::vector<double> > b(
            new std::vector<double>(std::move(
                common::import_vector_binary<double>(
                    params->basis_parameters()->basis_function_filename(a)
                    )
                    )
                    )); 
        return b;
    };
    import_wf = memoize(import_wf);

    std::function<PetscScalar (int, int)> findvalue = [&prototype,&params,&grid,&import_wf](int i, int j)->PetscScalar{
        if (i == j)
            return 0.0;
		auto a = import_wf(prototype[i]);
		//math::normalize(a,*grid);
		//std::cout << "a norm: " << 
		//std::cout << " after: " << math::normalize(a,*grid);
		auto b = import_wf(prototype[j]);
		//math::normalize(b,*grid);
		//std::cout << " b norm: " << 
		//std::cout << " after: " << math::normalize(a,*grid) << std::endl;
        PetscScalar radial = math::integrateGrid(
                *a,
                *b,
                *grid);
        PetscScalar angular = math::CGCoefficient<PetscScalar>(prototype[i],prototype[j]);
        //we are only considering spin up electrons.  so m_j always == +1/2 (since m_l is always 0)
        if (params->with_fs())
            angular *= std::sqrt ( 
                    (prototype[i].l + (prototype[i].j - 2 * prototype[i].l) * .5 * .5 + .5) * 
                    (prototype[j].l + (prototype[j].j - 2 * prototype[j].l) * .5 * .5 + .5) / 
                    ((2. * prototype[i].l + 1) * (2. * prototype[j].l + 1)));
       return radial * angular;
    };

   Mat H = common::populate_matrix(*params, 
                                    dipole_selection_rules, 
                                    findvalue, 
                                    prototype.size(),
                                    prototype.size());


   params->write_dipole_matrix(H);
   params->write_energy_eigenvalues();

   params->save_parameters();
   delete params;
   PetscFinalize();
   return 0;
}
