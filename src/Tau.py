try:
    from Recombination import *
except:
    from src.Recombination import *

def Find_fbh_Tau(
        Tau_max = 0.0561 + 2*0.0071,
        mc = 1e2,
        sbh = 1,
        zmin = 10,
        zmax = 500,
        nz = 2000,
        nm = 1000,
        Use_EoR = False,
        dmdm0_method = 0,
        Use_Halo_Boost = False,
        sbh_range = [-3, 5],
        show_status = False,
        Use_Halo_Boost_Interp = True,
        Use_Halo_Boost_Extrap = True,
        Use_mdot_Interp = True,
        Use_LowM_Extrap_mdot = True,
        Use_HighM_Extrap_mdot = True,
        Break_OmB_Lost_Ratio = 1.1,
        Fix_mdot = False,
        Fixed_mdot = 10,
        Show_Extrap_MSG = False,
        Ignore_Evo = False
    ):
    
    fbh_range = np.array([1e-30, 0.99])
    LgF_Range = np.log10(fbh_range)

    def dif_Tau(LgF):
        Rec = Call_HyRec(
            fbh0 = 10**LgF,
            mc = mc,
            sbh = sbh,
            Use_Ps_domain = False,
            zmin = zmin,
            zmax = zmax,
            nz = nz,
            nm = nm,
            Use_EoR = Use_EoR,
            dmdm0_method = dmdm0_method,
            Use_Halo_Boost = Use_Halo_Boost,
            sbh_range = sbh_range,
            show_status = False,
            Use_Halo_Boost_Interp = Use_Halo_Boost_Interp,
            Use_Halo_Boost_Extrap = Use_Halo_Boost_Extrap,
            Use_mdot_Interp = Use_mdot_Interp,
            Use_LowM_Extrap_mdot = Use_LowM_Extrap_mdot,
            Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot,
            Break_OmB_Lost_Ratio = Break_OmB_Lost_Ratio,
            Fix_mdot = Fix_mdot,
            Fixed_mdot = Fixed_mdot,
            Show_Extrap_MSG = Show_Extrap_MSG,
            Ignore_Evo = Ignore_Evo)
        zcut = 100
        z_vec = Rec['z'][::-1]
        xe_vec = Rec['xe'][::-1]
        tau = Find_Optical_Depth(z = z_vec, xe = xe_vec, xe_format = 0, zmax = zcut)
        r = tau - Tau_max
        return r
    try:
        LgF = Solve(F = dif_Tau, Xmin = LgF_Range[0], Xmax = LgF_Range[1], Precision = 1e-2, show_status = show_status)
        r = 10**LgF
    except:
        d0 = dif_Tau(LgF_Range[0])
        if d0 < 0:
            r = fbh_range[1]
        else:
            r = fbh_range[0]
    
    return r
