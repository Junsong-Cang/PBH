/* PBH accretion module
all masses are in msun unless specified otherwise
*/

// if Rhocr_C2_no_h2 is defined then this module is used by HyRec
#ifndef Rhocr_C2_no_h2
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "EFF_Tables.h"
#include "Useful_Functions.h"
#endif

// Ensure that the accretion rate does not exceed Eddintong limit
#define Allow_Super_Eddington 0
#define C_SI 299792458
#define Q_SI 1.602176634E-19
#define kB_SI 1.38064852E-23
#define mp_SI 1.67262158E-27
#define Solar_Mass 1.98847E30
#define Mbh_Evo_nz 100000
#define Radiation_Efficiency 0.1

double LCDM_Xe(double z)
{
  double r, f1, f2, zp;
  int id1, id2;
  zp = 1 + z;
  r = Interp_1D(zp, Redshift_Axis, LCDM_Xe_Template, Redshift_Size, 1, 1);
  return r;
}

double LCDM_Tm(double z)
{
  double r, f1, f2, zp;
  int id1, id2;
  zp = 1 + z;
  r = Interp_1D(zp, Redshift_Axis, LCDM_Tm_Template, Redshift_Size, 1, 1);
  return r;
}

double Hubble_LCDM(double z)
{
  double OmM, OmR, OmL, h, H0, r, zp;
  OmM = 0.3111;
  OmR = 9.065340263942517E-5;
  OmL = 1 - OmM - OmR;
  h = 0.6766;
  zp = 1 + z;
  H0 = 3.240755744239557e-18 * h;

  r = H0 * sqrt(OmL + OmM * pow(zp, 3) + OmR * pow(zp, 4));

  return r;
}

double Mdot_Eddington(double m)
{
  /* Eddington accretion rate, in kg/s
  ref : 2003.12589, Eq.10
  -- inputs --
  m : BH mass in msun
  */

  return 1.44E14 * m;
}

double Veff(double xe, double T, double z)
{
  // Effective speed in m/s
  // Eq.2&3 of 2108.11130
  double r, y, Cs, Vr;
  y = 5 / 3;
  Cs = sqrt(y * (1 + xe) * kB_SI * T / mp_SI);
  Vr = 30000 * fmin(1, (1 + z) / 1000);
  // printf("Cs = %E, Vr = %E\n", Cs, Vr);
  r = sqrt(pow(Cs, 2) + pow(Vr, 2));
  return r;
}

double viscosity(double m, double z, double xe, double T)
{
  // Gas viscosity withot DM accretion, Eq.3 of 2003.12589
  // m : bh mass in msun
  double r, zp, veff;
  zp = 1 + z;

  r = 0.257 + 1.45 * (xe / 0.01) * pow(zp / 1000, 2.5);
  r *= (m / 1E4) * pow(zp / 1000, 1.5);
  veff = Veff(xe, T, z);
  r *= pow(veff / 5740, -3);

  return r;
}

double Lambda(double m, double z, double xe, double T, double b)
{
  /* BH accretion efficiency
  ---- inputs ----
  b : viscosity
  */

  double r, xcr;
  xcr = (sqrt(1 + b) - 1) / b;
  r = xcr * xcr * exp(4.5 / (3 + pow(b, 0.75)));
  return r;
}

double K_factor(double m, double z, double xe, double T)
{
  // BH K factor, see 2003.12589
  double mh, zp, veff, r;
  zp = 1 + z;
  mh = 3 * m * 1000 / zp;
  veff = Veff(xe, T, z);
  r = 0.22 * (zp / 1000) * pow(mh, 2 / 3) * pow(veff / 1000, -2);

  return r;
}

double mdot_naked(double m, double z, double xe, double T)
{
  /* Dimensionless accretion rate for a nacked PBH
   */
  double r, l, vr, zp, b;
  b = viscosity(m, z, xe, T);
  vr = Veff(xe, T, z);
  l = Lambda(m, z, xe, T, b);
  zp = 1 + z;
  r = 0.4 * l * pow(zp / 1000, 3) * m * pow(vr / 1000, -3);

  return r;
}

double mdot_clothed(double m, double z, double xe, double T)
{
  /* Dimensionles accretion rate for clothed PBH
   */
  double k, r, zp, a, bh, b0, lh, y, veff, mdot, mh, p;
  a = 2.25; // halo density profile, 2108.11130

  zp = 1 + z;
  k = K_factor(m, z, xe, T);
  p = 2 - a;
  mh = 3 * m * 1000 / zp;
  if (k > 2)
  {
    r = mdot_naked(mh, z, xe, T);

    return r;
  }

  b0 = viscosity(m, z, xe, T);
  bh = b0 * pow(k, p / (1 - p));
  y = pow(1 + 10 * bh, 0.1);
  y *= exp(2 - k) * pow(k / 2, 2);
  lh = Lambda(m, z, xe, T, bh);
  lh *= pow(y, p / (1 - p));

  veff = Veff(xe, T, z);
  r = 0.023 * lh * (zp / 1000) * m * pow(veff / 5740, -3);

  return r;
}

double dMdz(double m, double z, double *DarkArray, int method, int model)
{
  /* Find dM/dz = -mdot/(1+z)/H
  unit : Msun
  -- inputs --
  *DarkArray : array containing info for xe and Tm
  method : how to get xe and T
           0 : Use LCDM template
           1 : Use updated DarkArray
  */
  double mdot, H, zp, xe, T, r;

  if (method == 0)
  {
    xe = LCDM_Xe(z);
    T = LCDM_Tm(z);
  }
  else
  {
    printf("Wrong method setting, DarkArray not supported yet.\n");
    exit(1);
  }
  if (model == 0)
  {
    mdot = mdot_clothed(m, z, xe, T);
  }
  else if (model == 1)
  {
    mdot = mdot_naked(m, z, xe, T);
  }
  H = Hubble_LCDM(z);
  zp = 1 + z;
  r = -mdot * Mdot_Eddington(m) / zp / H / Solar_Mass;
  return r;
}

double PBH_Accretion_Luminisoty(double m, double z, int model)
{
  /* Bomometric luminosity of an accreting PBH
  ignoring mass evolution for now
  unit : eV/s
  ---- inputs ----
  model : accretion model
          0 : naked
          1 : PBH clothed in DM halo
  */
  double xe, T, r, mdot, mdot_edd;
  // Don't use feedback for now
  xe = LCDM_Xe(z);
  T = LCDM_Tm(z);

  mdot = mdot_naked(m, z, xe, T);
  mdot_edd = Mdot_Eddington(m);

  r = mdot * mdot_edd * C_SI * C_SI / Q_SI;

  return r;
}

void Get_Mass_Evolution(double m0, double *z_out, double *m_out, int model)
{
  /* Find mass evolusion of a bh with initial mass m0 and store outputs in z_ou and m_out
  -- inputs --
  m0 : initial mass in msun
  z_out : output z axis, size = Mbh_Evo_nz
  m_out : output m axis, size = Mbh_Evo_nz
  */

  // internal settings
  double zmin = 0;
  double zmax = 3000;

  double vzp[Mbh_Evo_nz], m, dm, dz, z_prev, z;
  int id;
  Fill_Linear_Array(zmin + 1, zmax + 1, vzp, Mbh_Evo_nz, 1);

  m = m0;
  z_prev = zmax;

  for (id = Mbh_Evo_nz - 2; id > -1; id--)
  {
    z = vzp[id] - 1;
    dz = z - z_prev;
    dm = dMdz(m, z, vzp, 0, model) * dz;
    // dm = dMdz(m0, z, vzp, 0) * dz;
    m += dm;
    z_prev = z;
    printf("%E  %E  %E\n", z, m, m / m0);
  }
}

void Print_mdot(double m)
{
  double vz[1000];
  int nz = 100;
  int id;
  double mdot, z, xe, T;
  Fill_Linear_Array(5, 1000, vz, nz, 1);
  for (id = 0; id < nz; id++)
  {
    z = vz[id];
    xe = LCDM_Xe(z);
    T = LCDM_Tm(z);
    mdot = mdot_naked(m, z, xe, T);
    printf("%E  %E\n", z, mdot);
  }
}
