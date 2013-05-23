//ebss:
#include<common/common.hpp>
#include<common/math.hpp>
#include<common/parameters/MomentumParameters.hpp>
#include<common/special/coulomb/complex_functions.H>

//stl:
#include<vector>
#include<unordered_map>
#include<memory>
#include<functional>

//petsc:
#include<petsc.h>

template <typename ReturnType, typename Arg>
std::function<ReturnType (Arg)> memoize(std::function<ReturnType (Arg)> func)
{
    std::unordered_map<Arg, ReturnType> cache;
    return ([=](Arg t) mutable  {
            if (cache.find(t) == cache.end())
            {
                if (cache.size() > 200)
                    cache.clear();
                cache[t] = func(t);
            }
            return cache[t];
    });
}

int main(int argc, const char ** argv)
{
    using namespace std::placeholders;

    int ac = argc;
    char** av = new char*[argc];
    for (size_t i = 0; i < argc; i++)
    {
        av[i] = new char[strlen(argv[i])+1];
        std::copy(argv[i], argv[i] + strlen(argv[i])+1, av[i]);
    }
    PetscInitialize(&ac, &av, PETSC_NULL, PETSC_NULL);

    MomentumParameters<double> kparams(argc, argv, MPI_COMM_WORLD);
    if (kparams.rank() == 0)
    {
        std::cout << kparams.print();
    }

    auto hp = kparams.prototype();
    auto kp = kparams.kprototype();

    //for (auto a: kp)
        //std::cout << a << std::endl;
    //std::vector< std::vector<double> > wavefunctions(hp.size());
    //std::vector< std::vector<double> > expansion(kp.size());

    std::function< bool (int, int) > test = [&hp, &kp, &kparams](int i, int j) -> bool {
            return (hp[i].l == kp[j].l && hp[i].n < kparams.nmax() && hp[i].e.real() > 0.0); };

    std::function< std::shared_ptr<std::vector<double> > ( BasisID ) > import_wf = 
        [&kparams](BasisID a) -> std::shared_ptr<std::vector<double> > {
        std::shared_ptr< std::vector<double> > b(
            new std::vector<double>(std::move(
                common::import_vector_binary<double>(
                    kparams.hamiltonian().basis_parameters()->basis_function_filename(a)
                    )
                    )
                    )); 
        return b;
    };
    import_wf = memoize(import_wf);


    std::function< std::shared_ptr<std::vector<double> > (kBasisID ) > coulombf = 
        [&kparams](const kBasisID& a) -> std::shared_ptr<std::vector<double> > {
            try {
        std::shared_ptr< std::vector<double> > b(
                new std::vector<double>(std::move(
                    math::gsl_coulomb_wave_function( a, *(kparams.hamiltonian().basis_parameters()->grid()) )
                    )
                ));
        return b;
            } catch (std::out_of_range e) {
                std::shared_ptr< std::vector<double> > b(
                        new std::vector<double>(std::move(
                                math::coulomb_wave_function( a, *(kparams.hamiltonian().basis_parameters()->grid()) )
                                )
                            ));
        return b;
            }

    };
    coulombf = memoize(coulombf);

    std::function< std::complex<double> (int, int) > value = [&hp, &kp, &kparams, &import_wf, &coulombf](int i, int j){
            auto wf = import_wf(hp[i]);
            auto cwf = coulombf(kp[j]);
            return math::integrateTrapezoidRule(*wf, *cwf, *(kparams.hamiltonian().basis_parameters()->grid()), [](double r) -> double { return 1.; });
    };


    Mat etok = common::populate_matrix<std::complex<double> >( kparams, test, value, hp.size(), kp.size(), false);

    kparams.write_matrix( etok );
    kparams.save_parameters();

    PetscFinalize();
    return 0;
}
