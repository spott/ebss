#include <common/common.hpp>
#include <iostream>
#include <string>
#include <sstream>
#include <array>
#include <algorithm>

extern "C" {
#include <fftw3.h>
}
#include <complex>



int
main( int argc, const char ** argv )
{
    //takes six arguments: the min, the max, the dtau the base folder filename, and the output filename:
    if (argc < 6)
        return -1;

    double min, max, dtau;
    std::string rep_name;
    std::string output_fname;
    size_t N;

    min = stod( std::string( argv[1] ) );
    max = stod( std::string( argv[2] ) );
    dtau = stod( std::string( argv[3] ) );
    rep_name = std::string(argv[4]);
    output_fname = std::string(argv[5]);

    std::string first = rep_name.substr(0, rep_name.find("{"));
    std::string end = rep_name.substr(rep_name.find("}")+1);

    //std::cout << min << " " << max << " " << dtau << " " << rep_name << std::endl;
    std::cout << first << " " << end << std::endl;

    //import the files:
    std::vector< std::string > file_names;

    N = static_cast<size_t> ( ( max - min)/dtau);
    std::cout << "N: " << N << std::endl;

    char buffer[50];
    for (size_t i = 0; i < (max - min)/dtau; ++i)
    {
        sprintf( buffer, "%s%.1f%s/after.dat", first.c_str(), i * dtau, end.c_str() );
        //std::cout << buffer << "||" << std::endl;
        file_names.push_back(std::string( buffer ));
    }

    //for (auto & a : file_names )
        //std::cout << a << ";" << std::endl;

    //fftw stuff:
    double *in;
    fftw_complex *out;
    fftw_plan p;

    //malloc the stuff:
    in = (double*) fftw_malloc(sizeof(double) * N * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * (N/2 + 2));
    p = fftw_plan_dft_r2c_2d( N, N, in, out, FFTW_ESTIMATE);


    //move the data to *in.
    int current_place = 0;
    for (auto i = file_names.begin(); i < file_names.end(); ++i)
    {
        //import vector:
        std::cout << "importing file: " << *i << std::endl;
        auto v = common::import_vector_binary< std::array<double, 2> >(*i);

        auto dt = (v[1][0] - v[0][0]) * 2.41888432650516e-2;
        auto skip = static_cast<size_t> ( dtau / dt );

        if (skip * dt - dtau > 1e-10)
            std::cerr << "there isn't a good correlation between dtau and dt" << std::endl;
        std::cout << dt << " " << skip << std::endl;

        for (auto i = v.begin(); i - v.begin() < (skip * N); i += skip, ++current_place)
            in[current_place] = (*i)[1];

        std::cout << current_place << " " << std::endl;
    }
    std::cout << N*N << std::endl;

    std::cout << "starting the fft" << std::endl;
    fftw_execute(p);
    std::cout << "done with the fft" << std::endl;

    //save output:
    //create vector:
    std::vector< std::complex<double> > output( reinterpret_cast< std::complex<double>* >(out), reinterpret_cast< std::complex<double>* >(out + current_place ));
    common::export_vector_binary(output_fname, output);


    //cleanup
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
}
