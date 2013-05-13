//#include<common/new_parameters/BasisParameters.hpp>
//#include<common/new_parameters/HamiltonianParameters.hpp>
#include<common/parameters/PulsetrainParameters.hpp>
//#include<common/common.hpp>
#include<vector>
#include<array>
#include<petsc.h>

int
main (int argc, const char* argv[])
{
    int ac = argc;
    char** av = new char*[argc];
    for (size_t i = 0; i < argc; i++)
    {
        av[i] = new char[strlen(argv[i])];
        std::copy(argv[i], argv[i] + strlen(argv[i]), av[i]);
    }

    PetscInitialize(&ac,&av,PETSC_NULL,PETSC_NULL);

    //BasisParameters<double, double> p(argc, argv, PETSC_COMM_WORLD);
    //HamiltonianParameters h(argc, argv , PETSC_COMM_WORLD);
    PulsetrainParameters p(argc, argv, PETSC_COMM_WORLD);

    std::cout << p.print();

    double maxtime = p.max_time();
    std::cout << "maxtime: " << maxtime << std::endl;
    std::cout << "pulse_length: " << p.pulse_length() << std::endl;
    double t = 0;

    for (auto a: p.spacing())
        std::cout << a << ", ";
    std::cout << std::endl;

    std::vector< double > ef;
    std::vector< double > time;
    while (t < maxtime)
    {
        ef.push_back( p.efield(t).real() );
        time.push_back( t );
        //std::cerr << t << ", " << ef.back() << std::endl;
        t += p.dt();
    }
    common::export_vector_binary("efield.dat", ef);
    common::export_vector_binary("time.dat", time);

}
