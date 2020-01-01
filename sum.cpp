#include <cstdio>

#include "sumtk.hpp"

int main()
{
   const std::size_t dsize = 60000007;
   double* d = new double[dsize];

   for (std::size_t i = 0; i < dsize; ++i)
   {
      d[i] = (i & 1) ? 0.0000000123 : 12345.0000000001;
   }

   double sum = sumtk::sum(d,d + dsize,8);

   printf("sum: %35.15f\n",sum);

   delete[] d;
   return 0;
}
