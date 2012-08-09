#pragma once

#include<boost/function.hpp>
#include<boost/bind.hpp>


namespace numerov
{

   template<typename scalar>
   void find_basis(const parameters params,
                   const boost::function<scalar (scalar)> potential)
   {
      assert(params.type() == BASIS);

      PetscReal xmin = std::log(params.rmin());
      PetscReal xmax = std::log(params.rmax());
      PetscReal dx   = (xmax - xmin)/params.points();
      std::vector<PetscReal> grid(params.points());
      std::vector<PetscReal> pot(params.points());
      std::vector<PetscReal> wf(params.points());
      std::vector<PetscReal>::iterator git;
      std::vector<PetscReal>::iterator pit;
      std::vector<PetscReal>::iterator wfit;

      PetscReal energy = -.5;


      boost::function<scalar (scalar)> pot_x =
         boost::bind(potential, std::exp(_1));
      //setup grid and potential:
      for (git = grid.begin(), pit = pot.begin(); git < grid.end(); git++,pit++)
      {
         *git = xmin + (git - grid.begin()) * dx;
         *pit = std::exp(*git * 2) * (energy - potential(std::exp(*git)));
      }

      //bound!:  find split in the potential:
      std::vector<PetscReal>::iterator turnover = pot.begin();
      for (pit = pot.begin() + 1; pit < pot.end(); pit++)
      {
         if ((*(pit - 1) > 0 && *(pit) <= 0) ||
             (*(pit - 1) < 0 && *(pit) >= 0))
            turnover = pit;
         *pit = *pit * dx * dx / 12 + 1;
      }

      //choose first values
      *(wf.begin()) = 0;
      *(wf.begin()+1) = dx;
      *(wf.rend()) = 0;
      *(wf.rend()+1) = dx;
      //iterate using numerov's method
      numerov_step(wf.begin(), pot.begin(), turnover);
      numerov_step(wf.rend(), pot.rend(), turnover+1);

      //scale the meeting point:
      std::for_each(wf.begin(), turnover, )
   }

   template<typename scalar>
   void numerov_step(std::vector<scalar>::iterator wavefunction,
                     const std::vector<scalar>::iterator potential,
                     const std::vector<scalar>::iterator turnover)
   {
      //y_{n+1} = \frac {\left( 2-\frac{5 h^2}{6} f_n \right) y_n -
      //\left( 1+\frac{h^2}{12}f_{n-1}
      //\right)y_{n-1}}{1+\frac{h^2}{12}f_{n+1}}

      for (;potential <= turnover; potential++, wavefunction++)
      {
         *wavefunction = ((12 - 10 * (*(potential - 1)))*(*(wavefunction - 1)) -
                         *(potential - 2) * (*(wavefunction - 2)))/(*potential);
      }
   }

}
