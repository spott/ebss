#pragma once

#include<cassert>
#include<vector>
#include<boost/bind.hpp>
#include<boost/function.hpp>
#include<boost/ref.hpp>

#include<slepc.h>
#include<common/parameters.hpp>
#include<common/common.hpp>
//#include<common/special/bspline.hpp>


namespace finite_element
{

template<int order, typename element>
PetscReal find_value(
const parameters params,
   const boost::function<PetscReal (PetscReal)> potential,
   const std::vector<PetscReal> knots,
   const int i, const int j, const PetscReal dr = .01
   )
{
// Here I'm using BSplines, but theoretically, I could use any
// finitely bound spline.
element<order, PetscReal> bsi (&knots, i);
element<order, PetscReal> bsj (&knots, j);

PetscReal T, V, r, total;
total = 0;
for (int i = 0; i <= (params.rmax() - params.rmin())/dr; i++)
      {
         r = i * dr + params.rmin();
         // Kinetic Energy:
         T = bsi(r) * bsj.derivative(r, 2);
         // Potential Energy:
         V = bsi(r) * bsj(r) * potential(r);

         total += -.5 * T + V;
      }

      return total * dr;
   }

   template<int order>
   bool test(int i, int j)
   {
      return ( i - j <= order || j - i <= order );
   }


template <int order, typename element>
   void find_basis(const parameters params,
                        const boost::function<PetscReal (PetscReal)> potential)
   {
      //check to make sure the right params is being passed
      //assert(params.type() == BASIS);

      // Knot structure, we will start with a constant spacing:
      std::vector<PetscReal> knots(params.points());
      std::vector<PetscReal>::iterator it;
      // dr is the number of points - the number of knots at the beginning.
      PetscReal dr = (params.rmax() - params.rmin())
         /(params.points() - order * 2);

      for (it = knots.begin(); it < knots.end(); it++ )
      {
         if (it - knots.begin() <= order)
            *it = params.rmin();
         if (knots.end() - it <= order)
            *it = params.rmax();
         *it = params.rmin() + (it - knots.begin()) * dr;
      }
      boost::function< bool (int, int) > t = boost::cref(test<order>);
      boost::function< PetscReal (int, int) > fv =
         boost::bind<PetscReal>(find_value<order,element>,
                                boost::ref(params),
                                boost::cref(potential),
                                boost::cref(knots),
                                _1, _2, 0.01);
      Mat H = common::populate_matrix(params,
                                      t,
                                      fv,
                                      params.points() - order*2);

      Vec xr, xi;
      MatGetVecs(H, PETSC_NULL, &xr);
      MatGetVecs(H, PETSC_NULL, &xi);

      EPS eps;
      EPSCreate(params.comm(), &eps);

      EPSSetOperators(eps, H, PETSC_NULL);
      EPSSetProblemType(eps, EPS_HEP);
      EPSSetFromOptions(eps);
      EPSSolve(eps);

      PetscInt its, nconv,maxit, nev;
      PetscReal error, tol, re, im;
      PetscScalar kr, ki;
      EPSGetIterationNumber(eps,&its);
      PetscPrintf(params.comm()," Number of iterations of the method: %D\n",its);
      EPSGetDimensions(eps,&nev,PETSC_NULL,PETSC_NULL);
      PetscPrintf(params.comm()," Number of requested eigenvalues: %D\n",nev);
      EPSGetTolerances(eps,&tol,&maxit);
      PetscPrintf(params.comm()," Stopping condition: tol=%.4G, maxit=%D\n",tol,maxit);

      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Display solution and clean up
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      /*
         Get number of converged approximate eigenpairs
      */
      EPSGetConverged(eps,&nconv);
      PetscPrintf(params.comm()," Number of converged eigenpairs: %D\n\n",nconv);

      if (nconv>0) {
         /*
           Display eigenvalues and relative errors
         */
         PetscPrintf(params.comm()
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
               PetscPrintf(params.comm()," %9F%+9F j %12G\n",re,im,error);
            } else {
               PetscPrintf(params.comm(),"   %12F       %12G\n",re,error);
            }
         }
         PetscPrintf(params.comm(),"\n");


      }
   }


}
