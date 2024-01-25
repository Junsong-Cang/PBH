/*

A fast interface for python, to use this module do following

1. compile shared lib with either of these cmd
   gcc -shared -o mdot_interp.so mdot_interp.c
   icc -shared -o mdot_interp.so -fPIC mdot_interp.c
2. call from python like this::

import ctypes

main_path = '/Users/cangtao/cloud/GitHub/PBH/'

# Load the shared library
lib = main_path + 'src/c/mdot_interp.so'
c_function_lib = ctypes.CDLL(lib)

# Declare the function signature
Double = ctypes.c_double
Int = ctypes.c_int

# specify the name of function you want
c_function = c_function_lib.mdot_interp
# set input type
c_function.argtypes = (Double, Double, Double, Int, Int, Int, Int)
# set result type
c_function.restype = Double

# Call the C function
result = c_function(30.4, 497.6, 0.04897468161,1, 1, 1, 1)
print(result)

*/

#include <stdio.h>
#include <math.h>
#include "mdot_tab.h"

double Read_mdot(int zid, int mid, int EoR)
{
   double r;
   if (EoR == 0)
   {
      r = mdot_vec_LCDM[zid][mid];
   }
   else
   {
      r = mdot_vec_EoR[zid][mid];
   }
   return r;
}

int Find_Index(double x, double xmin, double xmax, int nx)
{
   // Find left index for x in an array log-distributed between [xmin, xmax]
   double dlx, lx1, lx2, lx;
   int idx;

   lx1 = log(xmin);
   lx2 = log(xmax);
   dlx = (lx2 - lx1) / ((double)nx - 1);
   lx = log(x);
   idx = (int)floor((lx - lx1) / dlx);

   // allow extrap
   if (idx < 0)
   {
      idx = 0;
   }
   if (idx > nx - 1)
   {
      idx = nx - 1;
   }

   return idx;
}

double mdot_interp(double z, double m, double OmB, int EoR, int Use_Edd_Limit, int Use_LowM_Extrap, int Use_HighM_Extrap)
{
   double x1, x2, f11, f12, f21, f22, F1, F2, F, x, Extrap_Redundancy, overflow, r, OmB_ratio;
   int zid1, zid2, mid1, mid2;

   // Even in the most extreme cases (mdot = 10, no OmB feedback), mass variation should not exceed 1e6
   Extrap_Redundancy = 1.0e6;

   zid1 = Find_Index(z, z_min, z_max, nz);
   mid1 = Find_Index(m, m_min, m_max, nm);
   zid2 = zid1 + 1;
   mid2 = mid1 + 1;

   f11 = Read_mdot(zid1, mid1, EoR);
   f12 = Read_mdot(zid1, mid2, EoR);
   f21 = Read_mdot(zid2, mid1, EoR);
   f22 = Read_mdot(zid2, mid2, EoR);

   // printf("zid1 = %d, zid2 = %d, mid1 = %d, mid2 = %d\n", zid1, zid2, mid1, mid2);
   // printf("f11 = %E, f12 = %E, f21 = %E, f22 = %E\n", f11, f12, f21, f22);

   if (ok_to_use_log_mdot == 1)
   {
      // Use log for mdot if possible, can improve precision. 
      // If some of mdot is close to 0 then log may give -inf which is not good. This is automatically checked by data/Tables_src/mdot_main.py
      
      f11 = log(f11);
      f12 = log(f12);
      f21 = log(f21);
      f22 = log(f22);
   }

   x1 = log(m_axis[mid1]);
   x2 = log(m_axis[mid2]);
   x = log(m);

   F1 = (f12 - f11) * (x - x1) / (x2 - x1) + f11;
   F2 = (f22 - f21) * (x - x1) / (x2 - x1) + f21;

   // now do z

   x1 = log(z_axis[zid1] + 1.0);
   x2 = log(z_axis[zid2] + 1.0);
   x = log(1.0 + z);
   F = (F2 - F1) * (x - x1) / (x2 - x1) + F1;

   if (ok_to_use_log_mdot == 1)
   {
      F = exp(F);
   }

   // decide whether to return F
   if (m < m_min)
   {
      overflow = m_min / m;
      if (Use_LowM_Extrap && (overflow < Extrap_Redundancy))
      {
         r = F;
      }
      else
      {
         r = -1.0; // a sign of things went wrong
      }
   }
   else if (m > m_max)
   {
      overflow = m / m_max;
      if (Use_HighM_Extrap && (overflow < Extrap_Redundancy))
      {
         r = F;
      }
      else
      {
         r = -1.0;
      }
   }
   else
   {
      r = F;
   }

   //printf("z1 = %E, z = %E, z2 = %E\n", z_axis[zid1], z, z_axis[zid2]);
   //printf("m1 = %E, m = %E, m2 = %E\n", m_axis[mid1], m, m_axis[mid2]);
   
   OmB_ratio = OmB/0.04897468161;
   r = OmB_ratio * r;

   if (Use_Edd_Limit && (r > 10.0))
   {
      r = 10.0;
   }

   return r;
}

/*
int main()
{
   double z, m, r;
   
   z = 30.4;
   m = 497.6;

   r = mdot_interp(z, m, 0.04897468161, 1, 1, 1, 1);
   printf("%E\n", r);

   return 0;
}
*/
