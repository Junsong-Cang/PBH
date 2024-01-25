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
#define C_SI 299792458           // Speed of light
#define Q_SI 1.602176634E-19     // Electron charge
#define kB_SI 1.38064852E-23     // Boltzmann Constant
#define mp_SI 1.67262158E-27     // Proton mass in kg
#define me_SI 9.1093837E-31      // Electron mass in kg
#define Sigma_Thomson 6.652E-29  // Thomson Cross Section in m^2
#define G_SI 6.6740831313131E-11 // Gravitational constant, in SI unit
#define Solar_Mass 1.98847E30    // Mass of sun in kg
#define Mbh_Evo_nz 100000

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

double Veff(double z, double xe, double T)
{
  // Effective speed in m/s
  // Eq.2&3 of 2108.11130
  double r, y, Cs, Vr;
  y = 5 / 3;
  Cs = sqrt(y * (1 + xe) * kB_SI * T / mp_SI);
  Vr = 30000 * fmin(1, (1 + z) / 1000);
  r = sqrt(pow(Cs, 2) + pow(Vr, 2));
  /*
  model used in 2303.06616
      if (Vr > Cs)
      {
        r = sqrt(Cs * Vr);
      }
      else
      {
        r = Cs;
      }
  */
  return r;
}

double Lambda(double m, double z, double xe, double T)
{
  double rbv, veff, rho_cmb, bc, yc, l_iso, l_ad, r;
  veff = Veff(z, xe, T);
  rbv = G_SI * m * Solar_Mass / pow(veff, 3);
  rho_cmb = 4.177995E-14 * pow(1 + z, 4); // cmb energy density in J/m^3
  l_iso = 1.120422267584516;
  l_ad = 0.11618950038622;
  bc = 4 * xe * Sigma_Thomson * rho_cmb * rbv / (3 * mp_SI * C_SI);
  yc = 8 * xe * Sigma_Thomson * rho_cmb * rbv / (3 * me_SI * C_SI * (1 + xe));

  // Step by step
  r = pow(yc, 2) / (88 + pow(yc, 2));
  r = l_ad + (l_iso - l_ad) * pow(r, 0.22);
  r = r / l_iso; // first line
  r *= exp(4.5 / (3 + pow(bc, 3 / 4)));
  r /= pow(sqrt(1 + bc) + 1, 2);

  return r;
}

double mdot_naked(double m, double z)
{
  /*BHL accretion rate in kg/s
   */
  double p4, l, rhob, obh2, veff, xe, T, m_si, r;
  p4 = 12.5663706143591; // 4*pi
  obh2 = 0.02242;
  xe = LCDM_Xe(z);
  T = LCDM_Tm(z);
  l = Lambda(m, z, xe, T);
  rhob = obh2 * 1.879E-26 * pow(1 + z, 3);
  veff = Veff(z, xe, T);
  m_si = m * Solar_Mass;
  r = p4 * l * rhob * pow(G_SI * m_si, 2) / pow(veff, 3);

  return r;
}

double PBH_Accretion_Luminisoty(double m, double z)
{
  /*
  PBH Accretion Luminisoty in eV/s
  */
  double mdot, eta, mdot_dimensionless, mdot_edd, mdot_radiated, r;

  mdot = mdot_naked(m, z);
  mdot_edd = Mdot_Eddington(m);
  mdot_dimensionless = mdot / mdot_edd;
  if (mdot_dimensionless > 0.1)
  {
    eta = 0.1;
  }
  else
  {
    eta = mdot_dimensionless;
  }

  mdot_radiated = eta * mdot;
  if (!Allow_Super_Eddington)
  {
    mdot_radiated = fmin(mdot_radiated, mdot_edd);
  }
  r = mdot_radiated * C_SI * C_SI / Q_SI;
  return r;
}

void Print_mdot(double m)
{
  // for testing
  double vz[1000];
  int nz = 100;
  int id;
  double mdot_edd, z, r;
  mdot_edd = Mdot_Eddington(m);

  Fill_Linear_Array(5, 1000, vz, nz, 1);
  for (id = 0; id < nz; id++)
  {
    z = vz[id];
    r = mdot_naked(m, z) / mdot_edd;
    printf("%E  %E\n", z, r);
  }
}

void Print_Luminosity_Array()
{
  /* Print luminosity array for matlab
  needed for deposition efficiency computation
  */
  // set array size
  int nm = 4000;
  double m1 = 1, m2 = 1E6, vm[4000];
  int zid, mid;
  double z, m, r;
  FILE *OutputFile, *M_Axis_File;
  remove("../DarkSide_src/Accreting_PBH/data/Luminosity_Table.txt");
  remove("../DarkSide_src/Accreting_PBH/data/Mbh_Axis.txt");
  OutputFile = fopen("../DarkSide_src/Accreting_PBH/data/Luminosity_Table.txt", "a");
  M_Axis_File = fopen("../DarkSide_src/Accreting_PBH/data/Mbh_Axis.txt", "a");
  Fill_Linear_Array(m1, m2, vm, nm, 1);
  for (zid = 0; zid < Redshift_Size; zid++)
  {
    for (mid = 0; mid < nm; mid++)
    {
      z = Redshift_Axis[zid] - 1;
      m = vm[mid];
      r = PBH_Accretion_Luminisoty(m, z);
      fprintf(OutputFile, "%E ", r);
    }
    fprintf(OutputFile, "\n");
  }
  for (mid = 0; mid < nm; mid++)
  {
    fprintf(M_Axis_File, "%E\n", vm[mid]);
  }
  fclose(OutputFile);
  fclose(M_Axis_File);
}
