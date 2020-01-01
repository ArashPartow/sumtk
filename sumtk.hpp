/*
 ********************************************************************
 *                    C++ Summation Toolkit                         *
 *                                                                  *
 * Author: Arash Partow (2002)                                      *
 * URL: http://www.partow.net/programming/sumtk/index.html          *
 *                                                                  *
 * Copyright notice:                                                *
 * Free use of the C++ Summation Toolkit library is permitted under *
 * the guidelines and in accordance with the most current version   *
 * of the Common Public License.                                    *
 * http://www.opensource.org/licenses/cpl1.0.php                    *
 *                                                                  *
 ********************************************************************
*/


#ifndef INCLUDE_SUMTK_HPP
#define INCLUDE_SUMTK_HPP


#include <algorithm>
#include <cstdio>
#include <cmath>
#include <deque>
#include <limits>
#include <vector>


namespace sumtk
{
   namespace details
   {
      template <typename T>
      inline T equal(const T v0, const T v1)
      {
         static const T epsilon = T(0.0000000001);
         return (std::abs(v0 - v1) <= (std::max(T(1),std::max(std::abs(v0),std::abs(v1))) * epsilon)) ? T(1) : T(0);
      }

      template <typename T>
      inline void kahan_sum(T& sum, T& error, T v)
      {
         T x = v - error;
         T y = sum + x;
         error = (y - sum) - x;
         sum = y;
      }

      template <typename T, typename InputIterator>
      inline T sum_block08(InputIterator begin, InputIterator end)
      {
         static const std::size_t block_size = 8;
         const std::size_t size = std::distance(begin,end);
         if (size >= block_size)
         {
            T total[8] = { T(0) };
            T error[8] = { T(0) };

            for (std::size_t i = 0; i < (size / block_size); ++i, begin += block_size)
            {
               const T v0 = *(begin + 0); const T v1 = *(begin + 1);
               const T v2 = *(begin + 2); const T v3 = *(begin + 3);
               const T v4 = *(begin + 4); const T v5 = *(begin + 5);
               const T v6 = *(begin + 6); const T v7 = *(begin + 7);
               kahan_sum(total[0],error[0],v0); kahan_sum(total[1],error[1],v1);
               kahan_sum(total[2],error[2],v2); kahan_sum(total[3],error[3],v3);
               kahan_sum(total[4],error[4],v4); kahan_sum(total[5],error[5],v5);
               kahan_sum(total[6],error[6],v6); kahan_sum(total[7],error[7],v7);
            }

            const std::size_t remainder = size % 8;
            if (remainder >= 1) kahan_sum(total[0],error[0],*(begin + 0));
            if (remainder >= 2) kahan_sum(total[1],error[1],*(begin + 1));
            if (remainder >= 3) kahan_sum(total[2],error[2],*(begin + 2));
            if (remainder >= 4) kahan_sum(total[3],error[3],*(begin + 3));
            if (remainder >= 5) kahan_sum(total[4],error[4],*(begin + 4));
            if (remainder >= 6) kahan_sum(total[5],error[5],*(begin + 5));
            if (remainder >= 7) kahan_sum(total[6],error[6],*(begin + 6));

            {
               T t = T(0);
               T e = T(0);
               kahan_sum(t,e,total[0]); kahan_sum(t,e,total[1]);
               kahan_sum(t,e,total[2]); kahan_sum(t,e,total[3]);
               kahan_sum(t,e,total[4]); kahan_sum(t,e,total[5]);
               kahan_sum(t,e,total[6]); kahan_sum(t,e,total[7]);
               return t;
            }
         }
         else
         {
            T total = T(0);
            T error = T(0);
            for (std::size_t i = 0; i < size; ++i, ++begin)
            {
               kahan_sum(total,error,*(begin));
            }
            return total;
         }
      }

      template <typename T, typename InputIterator>
      inline T sum_block16(InputIterator begin, InputIterator end)
      {
         static const std::size_t block_size = 16;
         const std::size_t size = std::distance(begin,end);
         if (size >= block_size)
         {
            T total[block_size] = { T(0) };
            T error[block_size] = { T(0) };

            for (std::size_t i = 0; i < (size / block_size); ++i, begin += block_size)
            {
               const T v00 = *(begin +  0); const T v01 = *(begin +  1);
               const T v02 = *(begin +  2); const T v03 = *(begin +  3);
               const T v04 = *(begin +  4); const T v05 = *(begin +  5);
               const T v06 = *(begin +  6); const T v07 = *(begin +  7);
               const T v08 = *(begin +  8); const T v09 = *(begin +  9);
               const T v10 = *(begin + 10); const T v11 = *(begin + 11);
               const T v12 = *(begin + 12); const T v13 = *(begin + 13);
               const T v14 = *(begin + 14); const T v15 = *(begin + 15);

               kahan_sum(total[ 0],error[ 0],v00); kahan_sum(total[ 1],error[ 1],v01);
               kahan_sum(total[ 2],error[ 2],v02); kahan_sum(total[ 3],error[ 3],v03);
               kahan_sum(total[ 4],error[ 4],v04); kahan_sum(total[ 5],error[ 5],v05);
               kahan_sum(total[ 6],error[ 6],v06); kahan_sum(total[ 7],error[ 7],v07);
               kahan_sum(total[ 8],error[ 8],v08); kahan_sum(total[ 9],error[ 9],v09);
               kahan_sum(total[10],error[10],v10); kahan_sum(total[11],error[11],v11);
               kahan_sum(total[12],error[12],v12); kahan_sum(total[13],error[13],v13);
               kahan_sum(total[14],error[14],v14); kahan_sum(total[15],error[15],v15);
            }

            const std::size_t remainder = size % block_size;
            if (remainder >=  1) kahan_sum(total[ 0],error[ 0],*(begin +  0));
            if (remainder >=  2) kahan_sum(total[ 1],error[ 1],*(begin +  1));
            if (remainder >=  3) kahan_sum(total[ 2],error[ 2],*(begin +  2));
            if (remainder >=  4) kahan_sum(total[ 3],error[ 3],*(begin +  3));
            if (remainder >=  5) kahan_sum(total[ 4],error[ 4],*(begin +  4));
            if (remainder >=  6) kahan_sum(total[ 5],error[ 5],*(begin +  5));
            if (remainder >=  7) kahan_sum(total[ 6],error[ 6],*(begin +  6));
            if (remainder >=  8) kahan_sum(total[ 7],error[ 7],*(begin +  7));
            if (remainder >=  9) kahan_sum(total[ 8],error[ 8],*(begin +  8));
            if (remainder >= 10) kahan_sum(total[ 9],error[ 9],*(begin +  9));
            if (remainder >= 11) kahan_sum(total[10],error[10],*(begin + 10));
            if (remainder >= 12) kahan_sum(total[11],error[11],*(begin + 11));
            if (remainder >= 13) kahan_sum(total[12],error[12],*(begin + 12));
            if (remainder >= 14) kahan_sum(total[13],error[13],*(begin + 13));
            if (remainder >= 15) kahan_sum(total[14],error[14],*(begin + 14));

            {
               T t = T(0);
               T e = T(0);
               kahan_sum(t,e,total[ 0]); kahan_sum(t,e,total[ 1]);
               kahan_sum(t,e,total[ 2]); kahan_sum(t,e,total[ 3]);
               kahan_sum(t,e,total[ 4]); kahan_sum(t,e,total[ 5]);
               kahan_sum(t,e,total[ 6]); kahan_sum(t,e,total[ 7]);
               kahan_sum(t,e,total[ 8]); kahan_sum(t,e,total[ 9]);
               kahan_sum(t,e,total[10]); kahan_sum(t,e,total[11]);
               kahan_sum(t,e,total[12]); kahan_sum(t,e,total[13]);
               kahan_sum(t,e,total[14]); kahan_sum(t,e,total[15]);
               return t;
            }
         }
         else
         {
            T total = T(0);
            T error = T(0);
            for (std::size_t i = 0; i < size; ++i, ++begin)
            {
               kahan_sum(total,error,*(begin));
            }
            return total;
         }
      }
   }

   inline double sum(const double* begin, const double* end, const std::size_t block_size = 8)
   {
      switch (block_size)
      {
         case  8 : return details::sum_block08<double>(begin,end);
         case 16 : return details::sum_block16<double>(begin,end);
         default : return std::numeric_limits<double>::quiet_NaN();
      }
   }

   inline float sum(const float* begin, const float* end, const std::size_t block_size = 8)
   {
      switch (block_size)
      {
         case  8 : return details::sum_block08<float>(begin,end);
         case 16 : return details::sum_block16<float>(begin,end);
         default : return std::numeric_limits<float>::quiet_NaN();
      }
   }

   template <typename Allocator>
   inline double sum(const std::vector<double,Allocator>& v, const std::size_t block_size = 8)
   {
      return sum(v.data(),v.data() + v.size(),block_size);
   }

   template <typename Allocator>
   inline float sum(const std::vector<float,Allocator>& v, const std::size_t block_size = 8)
   {
      return sum(v.data(),v.data() + v.size(),block_size);
   }

   template <typename Allocator>
   inline double sum(const std::deque<double,Allocator>& v, const std::size_t block_size = 8)
   {
      switch (block_size)
      {
         case  8 : return details::sum_block08<double>(v.begin(),v.end());
         case 16 : return details::sum_block16<double>(v.begin(),v.end());
         default : return std::numeric_limits<double>::quiet_NaN();
      }
   }

   template <typename Allocator>
   inline float sum(const std::deque<float,Allocator>& v, const std::size_t block_size = 8)
   {
      switch (block_size)
      {
         case  8 : return details::sum_block08<float>(v.begin(),v.end());
         case 16 : return details::sum_block16<float>(v.begin(),v.end());
         default : return std::numeric_limits<float>::quiet_NaN();
      }
   }

   inline bool run_test()
   {
      bool result = true;

      {
         const std::size_t dsize = 60000007;
         double* d = new double[dsize];

         for (std::size_t i = 0; i < dsize; ++i)
         {
            d[i] = (i & 1) ? 0.0000000123 : 12345.0000000001;
         }

         double d_sum8 = sum(d,d + dsize,8);

         if (!details::equal(d_sum8,370350049380.3720000373))
         {
            printf("[run_test01] - ERROR - Failed sum_block8. Expected:%30.10f\tResult:%30.10f\n",
                   370350049380.3720000373,
                   d_sum8);
            result = false;
         }

         double d_sum16 = sum(d,d + dsize,16);

         if (!details::equal(d_sum16,370350049380.3720000373))
         {
            printf("[run_test01] - ERROR - Failed sum_block16. Expected:%30.10f\tResult:%30.10f\n",
                   370350049380.3720000373,
                   d_sum16);
            result = false;
         }

         delete[] d;

         if (!result) return result;
      }

      {
         const std::size_t dsize = 60000007;
         std::vector<double> d(dsize);

         for (std::size_t i = 0; i < dsize; ++i)
         {
            d[i] = (i & 1) ? 0.0000000123 : 12345.0000000001;
         }

         double d_sum8 = sum(d,8);

         if (!details::equal(d_sum8,370350049380.3720000373))
         {
            printf("[run_test02] - ERROR - Failed sum_block8. Expected:%30.10f\tResult:%30.10f\n",
                   370350049380.3720000373,
                   d_sum8);
            result = false;
         }

         double d_sum16 = sum(d,16);

         if (!details::equal(d_sum16,370350049380.3720000373))
         {
            printf("[run_test02] - ERROR - Failed sum_block16. Expected:%30.10f\tResult:%30.10f\n",
                   370350049380.3720000373,
                   d_sum16);
            result = false;
         }

         if (!result) return result;
      }

      {
         const std::size_t dsize = 80000000;
         double* d = new double[dsize];

         for (std::size_t i = 0; i < dsize; ++i)
         {
            d[i] = (3.0 / 10000000000.0);
         }

         double d_sum8 = sum(d,d + dsize,8);

         if (!details::equal(d_sum8,0.024))
         {
            printf("[run_test03] - ERROR - Failed sum_block8. Expected:%30.10f\tResult:%30.10f\n",
                     0.024,
                     d_sum8);
            result = false;
         }

         double d_sum16 = sum(d,d + dsize,16);

         if (!details::equal(d_sum16,0.024))
         {
            printf("[run_test03] - ERROR - Failed sum_block16. Expected:%30.10f\tResult:%30.10f\n",
                     0.024,
                     d_sum16);
            result = false;
         }

         delete [] d;
         if (!result) return result;
      }

      {
         const std::size_t dsize = 80000000;
         double* d = new double[dsize];

         for (std::size_t i = 0; i < dsize; ++i)
         {
            if (i & 1)
               d[i] = (3.0 / 10000000000.0);
            else
               d[i] = 1234567890.123456789;
         }

         double d_sum8 = sum(d,d + dsize,8);

         if (!details::equal(d_sum8,49382715604938271.68))
         {
            printf("[run_test04] - ERROR - Failed sum_block8. Expected:%30.10f\tResult:%30.10f\n",
                     49382715604938271.68,
                     d_sum8);
            result = false;
         }

         double d_sum16 = sum(d,d + dsize,16);

         if (!details::equal(d_sum16,49382715604938271.68))
         {
            printf("[run_test04] - ERROR - Failed sum_block16. Expected:%30.10f\tResult:%30.10f\n",
                     49382715604938271.68,
                     d_sum16);
            result = false;
         }

         delete [] d;
         if (!result) return result;
      }

      {
         const std::size_t dsize = 80000000;
         double* d = new double[dsize];

         for (std::size_t i = 0; i < dsize; ++i)
         {
            if (i < (dsize/2))
               d[i] = +(1.0 / 100000000000.0);
            else
               d[i] = -(1.0 / 100000000000.0);
         }

         double d_sum8 = sum(d,d + dsize,8);

         if (!details::equal(d_sum8,0.0))
         {
            printf("[run_test05] - ERROR - Failed sum_block8. Expected:%30.10f\tResult:%30.10f\n",
                     0.0,
                     d_sum8);
            result = false;
         }

         double d_sum16 = sum(d,d + dsize,16);

         if (!details::equal(d_sum16,0.0))
         {
            printf("[run_test05] - ERROR - Failed sum_block16. Expected:%30.10f\tResult:%30.10f\n",
                     0.0,
                     d_sum16);
            result = false;
         }

         delete [] d;
         if (!result) return result;
      }

      {
         const std::size_t dsize = 10000000;
         double* d = new double[dsize];

         for (std::size_t i = 0; i < dsize; ++i)
         {
             d[i] = 0.123456789012345;
         }

         double d_sum8 = sum(d,d + dsize,8);

         if (!details::equal(d_sum8,1234567.89012345))
         {
            printf("[run_test06] - ERROR - Failed sum_block8. Expected:%30.10f\tResult:%30.10f\n",
                     1234567.89012345,
                     d_sum8);
            result = false;
         }

         double d_sum16 = sum(d,d + dsize,16);

         if (!details::equal(d_sum16,1234567.89012345))
         {
            printf("[run_test06] - ERROR - Failed sum_block16. Expected:%30.10f\tResult:%30.10f\n",
                     1234567.89012345,
                     d_sum16);
            result = false;
         }

         delete [] d;
         if (!result) return result;
      }

      // Future Test:
      //{
      //   double v[] = { 1.0, 1e100, 1.0, -1e100 };
      //   const std::size_t d_size = 10000 * (sizeof(v) / sizeof(double));
      //   double* d = new double[d_size];
      //
      //   for (std::size_t i = 0; i < d_size; ++i)
      //   {
      //      d[i] = v[i % (sizeof(v) / sizeof(double))];
      //   }
      //
      //   double d_sum8 = sum(d,d + d_size);
      //
      //   if (!equal(d_sum8,20000.0))
      //   {
      //      printf("[run_test??] - ERROR - Failed sum_block8. Expected:%30.10f\tResult:%30.10f\n",
      //               20000.0,
      //               d_sum8);
      //      result = false;
      //   }
      //}

      return true;
   }
}

#endif
