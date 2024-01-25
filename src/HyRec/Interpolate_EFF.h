// DM/PBH Deposition efficiency interpolation module

// include these files if compiled independently
// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
#include "EFF_Tables.h"
#include "Useful_Functions.h"

double Interp_EFF_DM_Annihilation(double Mdm, double z, int Dep_Channel, int Reaction_Channel)
{
  /* Interpolate deposition efficiencies for DM annihilation
  -- inputs --
  Mdm : DM mass in GeV
  z: redshfit
  Dep_Channel: Deposition Channel
          1 - HIon
          3 - LyA
          4 - Heat
  Reaction_Channel: Product that DM decay/annihilate into
          1 - Photon
          2 - Electron
          3 - Higgs
          4 - Muon
          5 - Tau
          6 - Q
          7 - CHarm
          8 - Bottom
          9 - Top
          10 - W
          11 - Z
          12 - Gluon
  */

  double Mdm_eV, x_in, EFF;
  Mdm_eV = Mdm * 1.0E9;
  // Validate inputs
  if ((Dep_Channel != 1) && (Dep_Channel != 3) && (Dep_Channel != 4))
  {
    printf("Error: wrong choice of deposition channels.\n");
    exit(1);
  }

  // ---- Gamma ----
  if (Reaction_Channel == 1)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Phot_HIon_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Phot_LyA_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Phot_Heat_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Elec ----
  else if (Reaction_Channel == 2)
  {
    x_in = Mdm_eV - Electron_Mass_eV;
    if (x_in < 0.0)
    {
      printf("Error: DM mass too small to annihilate into electron+positron.\n");
      exit(1);
    }
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Elec_HIon_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Elec_LyA_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Elec_Heat_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Higgs ----
  else if (Reaction_Channel == 3)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Higgs_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Higgs_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Higgs_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Muon ----
  else if (Reaction_Channel == 4)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Muon_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Muon_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Muon_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Tau ----
  else if (Reaction_Channel == 5)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Tau_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Tau_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Tau_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Q ----
  else if (Reaction_Channel == 6)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Q_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Q_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Q_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Charm ----
  else if (Reaction_Channel == 7)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Charm_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Charm_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Charm_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Bottom ----
  else if (Reaction_Channel == 8)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Bottom_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Bottom_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Bottom_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Top ----
  else if (Reaction_Channel == 9)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Top_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Top_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Top_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- W ----
  else if (Reaction_Channel == 10)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_W_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_W_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_W_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Z ----
  else if (Reaction_Channel == 11)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Z_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Z_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Z_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Gluon ----
  else if (Reaction_Channel == 12)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Gluon_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Gluon_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Gluon_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  else
  {
    printf("Unknown DM decay/annihilation channel.\n");
    exit(1);
  }
  return EFF;
}

double Interp_EFF_DM_Decay(double Mdm, double z, int Dep_Channel, int Reaction_Channel)
{
  /* Interpolate deposition efficiencies for DM decay
  -- inputs --
  Mdm : DM mass in GeV
  z: redshfit
  Dep_Channel: Deposition Channel
          1 - HIon
          3 - LyA
          4 - Heat
  Reaction_Channel: Product that DM decay/annihilate into
          1 - Photon
          2 - Electron
          3 - Higgs
          4 - Muon
          5 - Tau
          6 - Q
          7 - CHarm
          8 - Bottom
          9 - Top
          10 - W
          11 - Z
          12 - Gluon
  */

  double Mdm_eV, x_in, EFF;
  Mdm_eV = Mdm * 1.0E9;
  // Validate inputs
  if ((Dep_Channel != 1) && (Dep_Channel != 3) && (Dep_Channel != 4))
  {
    printf("Error: wrong choice of deposition channels\n");
    exit(1);
  }

  // ---- Gamma ----
  if (Reaction_Channel == 1)
  {
    x_in = Mdm_eV / 2.0;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Phot_HIon_Decay, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Phot_LyA_Decay, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Phot_Heat_Decay, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Electron ----
  else if (Reaction_Channel == 2)
  {
    x_in = Mdm_eV / 2.0 - Electron_Mass_eV;
    if (x_in < 0.0)
    {
      printf("Error: DM mass too small to decay into electron+positron.\n");
      exit(1);
    }

    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Elec_HIon_Decay, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Elec_LyA_Decay, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Elec_Heat_Decay, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Higgs ----
  else if (Reaction_Channel == 3)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Higgs_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Higgs_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Higgs_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Muon ----
  else if (Reaction_Channel == 4)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Muon_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Muon_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Muon_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Tau ----
  else if (Reaction_Channel == 5)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Tau_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Tau_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Tau_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Q ----
  else if (Reaction_Channel == 6)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Q_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Q_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Q_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Charm ----
  else if (Reaction_Channel == 7)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Charm_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Charm_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Charm_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Bottom ----
  else if (Reaction_Channel == 8)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Bottom_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Bottom_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Bottom_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Top ----
  else if (Reaction_Channel == 9)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Top_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Top_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Top_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- W ----
  else if (Reaction_Channel == 10)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_W_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_W_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_W_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Z ----
  else if (Reaction_Channel == 11)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Z_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Z_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Z_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Gluon ----
  else if (Reaction_Channel == 12)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Gluon_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Gluon_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Gluon_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  else
  {
    printf("Unknown DM decay/annihilation channel.\n");
    exit(1);
  }

  return EFF;
}

int Convert_to_Int(double x)
{
  int ix;
  ix = (int)round(x);
  return ix;
}

double Interp_EFF_Hawking(double Mbh, double z, double PBH_Spin, double PBH_Model, int Dep_Channel)
{
  /* Interpolate PBH Hawking Radiation deposition efficiencies
  inputs:
  Mbh : PBH mass in gram
  z : redshift
  PBH_Spin : PBH_Spin
  PBH_Model : PBH_Model defined in energy_injection.h
  Dep_Channel: Deposition Channel
          1 - HIon
          3 - LyA
          4 - Heat
  */
  double EFF, Tbh;
  int Spin_id, model_id;
  model_id = Convert_to_Int(PBH_Model);
  Tbh = 1.06E6 * (1.0E16 / Mbh); // BH temperature in eV
  if (PBH_Spin < 0.2)
  {
    Spin_id = 1;
  }
  else if ((PBH_Spin > 0.24) && (PBH_Spin < 0.26))
  {
    Spin_id = 2;
  }
  else if ((PBH_Spin > 0.49) && (PBH_Spin < 0.51))
  {
    Spin_id = 3;
  }
  else if ((PBH_Spin > 0.74) && (PBH_Spin < 0.76))
  {
    Spin_id = 4;
  }
  else if ((PBH_Spin > 0.991) && (PBH_Spin < 0.9995))
  {
    Spin_id = 5;
  }
  else if ((PBH_Spin > 0.9996) && (PBH_Spin < 0.99995))
  {
    Spin_id = 6;
  }

  if (Spin_id == 1)
  {
    if (model_id == 2)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K1, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K1, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K1, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
    else if (model_id == 3)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K1B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K1B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K1B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
  }
  else if (Spin_id == 2)
  {
    if (model_id == 2)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K2, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K2, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K2, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
    else if (model_id == 3)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K2B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K2B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K2B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
  }
  else if (Spin_id == 3)
  {
    if (model_id == 2)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K3, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K3, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K3, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
    else if (model_id == 3)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K3B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K3B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K3B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
  }
  else if (Spin_id == 4)
  {
    if (model_id == 2)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K4, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K4, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K4, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
    else if (model_id == 3)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K4B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K4B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K4B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
  }
  else if (Spin_id == 5)
  {
    if (model_id == 2)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K5, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K5, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K5, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
    else if (model_id == 3)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K5B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K5B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K5B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
  }
  else if (Spin_id == 6)
  {
    if (model_id == 2)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K6, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K6, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K6, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
    else if (model_id == 3)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K6B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K6B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K6B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
  }

  return EFF;
}

double Interp_EFF_Accreting_PBH(double m, double z, int Dep_Channel)
{
  /* Interpolate deposition efficiencies for accreting PBH
  -- inputs --
  m : PBH mass in msun
  */

  double r;
  // Validate inputs
  if ((Dep_Channel != 1) && (Dep_Channel != 3) && (Dep_Channel != 4))
  {
    printf("Error: wrong choice of deposition channels.\n");
    exit(1);
  }
  if (Dep_Channel == 1)
  {
    r = Interp_2D(m, 1 + z, PBH_Accretion_Mass_Axis, Redshift_Axis, EFF_Accreting_PBH_Naked_HIon, PBH_Accretion_Mass_Axis_Size, Redshift_Size, 1);
  }
  else if (Dep_Channel == 3)
  {
    r = Interp_2D(m, 1 + z, PBH_Accretion_Mass_Axis, Redshift_Axis, EFF_Accreting_PBH_Naked_LyA, PBH_Accretion_Mass_Axis_Size, Redshift_Size, 1);
  }
  else if (Dep_Channel == 4)
  {
    r = Interp_2D(m, 1 + z, PBH_Accretion_Mass_Axis, Redshift_Axis, EFF_Accreting_PBH_Naked_Heat, PBH_Accretion_Mass_Axis_Size, Redshift_Size, 1);
  }
  return r;
}
