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

struct rid {
    double r;
    int l;
}


int main(int argc, const char ** argv)
{
    int ac = argc;
    char** av = new char*[argc];
    for (size_t i = 0; i < argc; i++)
    {
        av[i] = new char[strlen(argv[i])+1];
        std::copy(argv[i], argv[i] + strlen(argv[i])+1, av[i]);
    }
    PetscInitialize(&ac, &av, PETSC_NULL, PETSC_NULL);

    int npoints = 10000;
    int lmax = 30;
    double gridmax = 1000.;

    std::vector<rid> rprototype(npoints * lmax);
    for ( auto i = rprototype.begin(); i < rprototype.end(); ++i)
    {
        for 

    std::function< PetscScalar (int, int) > test = [](int i, int j)

}
