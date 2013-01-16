#pragma once

template <typename write_type_ = double>
class MomentumParameters: public Parameters
{
public:
    MomentumParameters(int argc, const char** argv, MPI_Comm comm): Parameters(comm)
    {
        register_parameters();

        opt.parse(argc, argv);

        if (opt.isSet("+d")) {
            std::string pretty;
            opt.prettyPrint(pretty);
            std::cout << pretty;
        }
        if (opt.isSet("-h")) {
            std::string usage;
            opt.getUsage(usage,80,ez::ezOptionParser::ALIGN);
            std::cout << usage;
        }

        if (opt.isSet("-momentum_config"))
        {
            std::string fname;
            opt.get("-hamiltonian_config")->getString(fname);
            if (! opt.importFile(fname.c_str(), '#'))
            {
                std::cout << "file must exist!" << std::endl;
                throw std::exception();
            }
        }

        opt.get("-hamiltonian_nmax")->getInt(nmax_);
        opt.get("-hamiltonian_lmax")->getInt(lmax_);
        opt.get("-hamiltonian_mmax")->getInt(mmax_);
        opt.get("-hamiltonian_folder")->getString(folder_);


        folder_ = common::absolute_path(folder_);

