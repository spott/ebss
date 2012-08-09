#pragma once

#include<vector>


namespace BSpline
{
template <int order, typename scalar>
class BSpline
{
public:
   BSpline(const std::vector<scalar> *knots_,
           const int num_, const scalar scalingFactor_=1)
      : knots(knots_), num(num_), scalingFactor(scalingFactor_)
   {}

   int left_knot(scalar x)
   {
      for (int i = 1; i < (*knots).size(); i++)
         if ((*knots)[i-1] <= x && x < (*knots)[i])
            return i-1;
      return (*knots).back();
   }

   scalar derivative(scalar x, int dn)
   {
      scalar saved = 0;
      scalar temp = 0;
      scalar Uleft = 0;
      scalar Uright = 0;
      scalar deriv[dn+1];
      scalar N[order+1][order+1];
      scalar ND[order+1];
      std::vector<scalar> U = *knots;
      if (x < U[num] || x >= U[num+order+1])
         return 0;
      for (int j = 0; j <= order; j++)
      {
         if ( x >= U[num+j] && x < U[num + j + 1])
            N[j][0] = 1;
         else
            N[j][0] = 0;
      }
      for (int k = 1; k <= order; k++)
      {
         if (N[0][k-1] == 0)
            saved = N[0][k-1];
         else
            saved = ((x-U[num]) * N[0][k-1]) / (U[num + k] - U[num]);
         for (int j = 0; j< order - k + 1; j++)
         {
            Uleft = U[num + j + 1];
            Uright = U[num + j + k + 1];
            if (N[j+1][k-1] == 0)
            {
               N[j][k] = saved;
               saved = 0;
            }
            else
            {
               temp = N[j+1][k-1]/(Uright-Uleft);
               N[j][k] = saved + (Uright - x) * temp;
               saved = (x - Uleft) * temp;
            }
         }
      }
      deriv[0] = N[0][order];
      for (int k = 1; k <= dn; k++)
      {
         for (int j = 0; j<= k; j++)
            ND[j] = N[j][order-k];
         for (int jj = 1; jj <= k; jj++)
         {
            if (ND[0] == 0)
               saved = 0;
            else
               saved = ND[0] / (U[num + order - k + jj] - U[num]);
            for (int j = 0; j < k - jj + 1; j++)
            {
               Uleft = U[num+j+1];
               Uright = U[num + j + order + jj + 1];
               if (ND[j+1] == 0)
               {
                  ND[j] = (order - k + jj) * saved;
                  saved = 0;
               }
               else
               {
                  temp = ND[j+1]/(Uright - Uleft);
                  ND[j] = (order-k+jj)*(saved - temp);
                  saved = temp;
               }
            }
         }
         deriv[k] = ND[0];
      }
      return deriv[dn];
   }

   scalar operator()(scalar x)
   {
      scalar saved = 0;
      scalar temp = 0;
      scalar Uleft = 0;
      scalar Uright = 0;
      std::vector<scalar> N(order + 1);
      if (( num == 0 && x == (*knots)[0]) ||
          (x == knots->size() - order - 1 && x == knots->back()))
         return 1.0;
      if ( x < (*knots)[num] || x >= (*knots)[num + order + 1])
         return 0;
      for (int j = 0; j <= order; j++)
      {
         if (x >= (*knots)[num + j] && x < (*knots)[num + j + 1])
            N[j] = 1.0;
         else
            N[j] = 0;
      }
      for (int k = 1; k <= order; k++)
      {
         if (N[0] == 0)
            saved = 0;
         else
            saved = ((x - (*knots)[num])*N[0]) /
               ((*knots)[num+k] - (*knots)[num]);
         for (int j = 0; j < order - k + 1; j++)
         {
            Uleft = (*knots)[num + j + 1];
            Uright = (*knots)[num + j + k + 1];
            if (N[j+1] == 0)
            {
               N[j] = saved;
               saved = 0;
            }
            else
            {
               temp = N[j+1]/(Uright-Uleft);
               N[j] = saved + (Uright - x)*temp;
               saved = (x - Uleft) * temp;
            }
         }
      }
      return N[0] * scalingFactor;
   }

private:
   const std::vector<scalar> *knots;
   const int num;
   const scalar scalingFactor;
};

}
