#pragma once

#include<cassert>
#include<vector>

#include<slepc.h>
#include<common/parameters/Parameters.hpp>
#include<common/parameters/BasisParameters.hpp>
#include<common/common.hpp>


namespace finite_difference
{
    //template<int order, typename scalar>

        //std::vector<BasisID> find_basis_l( std::function<scalar (scalar)> potential,
                                           //BasisParameters<scalar>* params, 
                                           //std::vector<scalar>* grid)
        //{





        template<int order, typename scalar>
        void find_basis(std::function<scalar (scalar)> pot,
                BasisParameters<scalar>* params)
        {
            // dr is the number of points - the number of grid at the beginning.
            scalar dr = (params->rmax() - params->rmin())
                /(params->points());

            std::cout << order << std::endl;

            std::vector<scalar> *rgrid = params->grid();
            for (size_t i = 0; i < rgrid->size(); i++)
                rgrid->at(i) =params->rmin() + i * dr;
            std::cout << rgrid->size() << std::endl;

            std::vector<BasisID> *energies = params->basis_prototype();
            energies->resize(0);

            std::function< bool (int, int) > t = [](int i, int j) { return ( std::abs(i - j) < order || std::abs(j - i) < order ); };
            BasisID tmp;
            for (int l = 0; l < params->lmax(); l++)
            {
                std::function< scalar ( scalar ) > potential = [pot,l](scalar r) { return pot(r) + l * (l+1) / (2 * r * r); };

                std::function< scalar (int, int) > fv = [potential, dr, rgrid] (int i, int j)
                {
                    if (i == j) // diagonal
                        return potential(rgrid->at(i)) + 1 / (dr*dr);
                    else if ( i - 1 == j || i + 1 == j) //off diagonal
                        return -1/(2 * dr*dr);
                };

                Mat H = common::populate_matrix( static_cast<const Parameters>(*params),
                        t,
                        fv,
                        params->points(),
                        params->points(),
                        true);

                Vec xr, xi;
                MatGetVecs(H, PETSC_NULL, &xr);
                MatGetVecs(H, PETSC_NULL, &xi);

                EPS eps;
                EPSCreate(params->comm(), &eps);

                EPSSetOperators(eps, H, PETSC_NULL);
                EPSSetProblemType(eps, EPS_HEP);
                EPSSetDimensions(eps, params->nmax() - l, PETSC_DECIDE, PETSC_DECIDE);
                EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL);
                EPSSetFromOptions(eps);
                if (params->rank() == 0) std::cout << "starting solve" << std::endl;
                EPSSolve(eps);

                PetscInt its, nconv,maxit, nev;
                PetscReal error, tol, re, im;
                PetscScalar kr, ki;
                EPSGetIterationNumber(eps,&its);
                PetscPrintf(params->comm()," Number of iterations of the method: %D\n",its);
                EPSGetDimensions(eps,&nev,PETSC_NULL,PETSC_NULL);
                PetscPrintf(params->comm()," Number of requested eigenvalues: %D\n",nev);
                EPSGetTolerances(eps,&tol,&maxit);
                PetscPrintf(params->comm()," Stopping condition: tol=%.4G, maxit=%D\n",tol,maxit);

                /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                   Display solution and clean up
                   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
                /*
                   Get number of converged approximate eigenpairs
                   */
                EPSGetConverged(eps,&nconv);
                PetscPrintf(params->comm()," Number of converged eigenpairs: %D\n\n",nconv);

                if (nconv>0) {
                    /*
                       Display eigenvalues and relative errors
                       */
                    PetscPrintf(params->comm()
                            ,
                            "           k          ||Ax-kx||/||kx||\n"
                            "   ----------------- ------------------\n");

                    for (int i=0;i<nconv;i++) {
                        /*
                           Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
                           ki (imaginary part)
                           */
                        EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);
                        /*
                           Compute the relative error associated to each eigenpair
                           */
                        EPSComputeRelativeError(eps,i,&error);

                        re = PetscRealPart(kr);
                        im = PetscImaginaryPart(kr);



                        if (im!=0.0) {
                            PetscPrintf(params->comm()," %9F%+9F j %12G\n",re,im,error);
                        } else {
                            PetscPrintf(params->comm(),"   %12F       %12G\n",re,error);
                        }
                        energies->push_back({i+1+l,l,0,kr});
                        auto vout1 = common::Vec_to_vector(xr);
                        if (params->rank() == 0)
                            common::export_vector_binary<PetscReal>(
                                    params->basis_function_filename({i+1+l,l,0,kr}), 
                                    common::vector_type_change<PetscScalar, PetscReal>(vout1));
                    }
                    PetscPrintf(params->comm(),"\n");
                }
            

                MatDestroy(&H);
                EPSDestroy(&eps);
                VecDestroy(&xr);
                VecDestroy(&xi);
            }
        }


}
