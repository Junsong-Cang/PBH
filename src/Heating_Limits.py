try:
    from Recombination import *
except:
    from src.Recombination import *

def Find_fbh_Tk(
        mc = 1e2,
        sbh = 1,
        z_for_tk_max = [4.3, 4.8, 6.08],
        tk_max = [20000, 10000.0, 18620.87136],
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

    def Tk_dif(LgF):
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
        
        z_vec = Rec['z'][::-1]
        t_vec = Rec['t'][::-1]
        z_axis = np.array(z_for_tk_max)
        t_max = np.array(tk_max)
        t_model = np.interp(x = z_axis, xp = z_vec, fp = t_vec)
        dif = t_model > t_max
        if True in dif:
            r = 1.0
        else:
            r = -1.0
        return r
    
    try:
        LgF = Solve(F = Tk_dif, Xmin = LgF_Range[0], Xmax = LgF_Range[1], Precision = 1e-2, show_status = show_status)
        r = 10**LgF
    except:
        d0 = Tk_dif(LgF_Range[0])
        d1 = Tk_dif(LgF_Range[1])
        if d0 < 0:
            r = fbh_range[1]
        else:
            r = fbh_range[0]
    return r
