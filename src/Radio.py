# Radio excess module

try:
    from src.main import *
except:
    from main import *

def Comoving_Radio_Emissivity(
        z_axis, 
        EMS_axis,
        z = np.logspace(1,2,100),
        freq = 1.0):
    '''
    Find comoving radio emissivity from given EMS axis
    ----inputs----
    z : redshift
    freq : frequency in GHz
    EMS_axis : emissivity data @ 1GHz
    z_axis : z axis for EMS_axis
    '''
    aR = 0.62

    zp = 1+z

    # normaly z_axis should be decreasing but let's check just in case
    if z_axis[0] > z_axis[-1]:
        zp_vec = z_axis[::-1] + 1
        EMS_vec = EMS_axis[::-1]
    else:
        zp_vec = z_axis + 1
        EMS_vec = EMS_axis

    EMS = np.interp(x = np.log(zp), xp = np.log(zp_vec), fp = EMS_vec, left = 0.0, right = 0.0) # for 1 GHz
    r = EMS * freq**(-aR)

    return r

def T_Radio_fun(
        z_axis,
        EMS_axis,
        z = 10):
    '''
    Compute radio brightness temp @ 21cm frequency
    not very efficient because H is computted too often but I guess this is very fast anyway
    '''
    nz = 1000 # should be enough
    Prefix = 3.8025476052650717e+28 # c^3/(8 * pi * kB * v21^2), in SI unit
    zp3 = (1+z)**3
    z1 = z
    z2 = np.max(z_axis)
    z_prev = np.logspace(np.log10(z1 + 1), np.log10(z2 + 1), nz) - 1
    v21 = 1.429 # GHz
    v = v21 * (1+z_prev)/(1+z)
    EMS = Comoving_Radio_Emissivity(z_axis = z_axis, EMS_axis = EMS_axis,z = z_prev, freq = v)
    H = Hubble(z_prev)
    
    fun = zp3 * Prefix * EMS / ((1+z_prev) * H)

    r = np.trapz(x = z_prev, y = fun)

    return r

def T_Radio_Kernel(z_axis, EMS_axis):
    '''
    Get T_Radio from 0 to max(z_axis)
    '''
    nz = 1000
    z1 = 0
    z2 = np.max(z_axis)
    z_vec = np.logspace(np.log10(z1 + 1), np.log10(z2 + 1), nz) - 1
    T = np.zeros(nz)
    
    for idx in np.arange(0, nz):
        z = z_vec[idx]
        T[idx] = T_Radio_fun(z_axis = z_axis, EMS_axis = EMS_axis, z = z)
    
    r = {'z' : z_vec,
        'T' : T} # T_Radio @ 21cm
    
    return r
    
def T_Radio(
        mc = 1e2,
        sbh = 1,
        fbh0 = 1e-5,
        Use_Ps_domain = False,
        PsA = 0.1,
        kc = 1e5,
        SigmaPs = 0.5,
        DeltaC = 0.45,
        OmB0 = cosmo['OmB'],
        Use_EoR = True,
        dmdm0_method = 0,
        Use_Halo_Boost = False,
        zmin = 8,
        zmax = 1000,
        nz = 1500,
        nm = 500,
        sbh_range = [-3, 5],
        show_status = False,
        Use_Halo_Boost_Interp = True,
        Use_Halo_Boost_Extrap = True,
        Use_mdot_Interp = True,
        Use_LowM_Extrap_mdot = True,
        Use_HighM_Extrap_mdot = True,
        Break_OmB_Lost_Ratio = 1.1,
        Show_Extrap_MSG = False,
        Fix_mdot = False,
        Fixed_mdot = 10,
        Ignore_Evo = False
    ):
    
    if sbh < 0:
        Evo = Get_Evo_Lite(
            fbh0 = fbh0,
            m0 = mc,
            OmB0 = OmB0,
            zmin = zmin,
            zmax = zmax,
            nz = nz,
            Use_Halo_Boost = Use_Halo_Boost,
            Use_EoR = Use_EoR,
            Use_OmB_FeedBack = True,
            Fix_mdot = Fix_mdot,
            Fixed_mdot = Fixed_mdot,
            Use_Halo_Boost_Interp = Use_Halo_Boost_Interp,
            Use_Halo_Boost_Extrap = Use_Halo_Boost_Extrap,
            Use_mdot_Interp = Use_mdot_Interp,
            Use_LowM_Extrap_mdot = Use_LowM_Extrap_mdot,
            Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot,
            Show_Extrap_MSG = Show_Extrap_MSG,
            Ignore_Evo = Ignore_Evo)
        z_axis = Evo['z']
        EMS_axis = Evo['Radio_EMS']
    else:
        Evo =  Get_Evo(
            mc = mc,
            sbh = sbh,
            fbh0 = fbh0,
            Use_Ps_domain = Use_Ps_domain,
            PsA = PsA,
            kc = kc,
            SigmaPs = SigmaPs,
            DeltaC = DeltaC,
            OmB0 = OmB0,
            Use_EoR = Use_EoR,
            dmdm0_method = dmdm0_method,
            Use_Halo_Boost = Use_Halo_Boost,
            zmin = zmin,
            zmax = zmax,
            Use_HZ = False,
            nz = nz,
            nm = nm,
            sbh_range = sbh_range,
            show_status = show_status,
            Use_Halo_Boost_Interp = Use_Halo_Boost_Interp,
            Use_Halo_Boost_Extrap = Use_Halo_Boost_Extrap,
            Use_mdot_Interp = Use_mdot_Interp,
            Use_LowM_Extrap_mdot = Use_LowM_Extrap_mdot,
            Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot,
            Break_OmB_Lost_Ratio = Break_OmB_Lost_Ratio,
            Show_Extrap_MSG = Show_Extrap_MSG,
            Fix_mdot = Fix_mdot,
            Fixed_mdot = Fixed_mdot,
            Ignore_Evo = Ignore_Evo)
        z_axis = Evo['z_vec']
        EMS_axis = Evo['Radio_EMS']

    r = T_Radio_Kernel(z_axis = z_axis, EMS_axis = EMS_axis)
    
    return r

def Find_fbh_Radio(
        mc = 1e2,
        sbh = 1,
        OmB0 = cosmo['OmB'],
        Use_EoR = True,
        dmdm0_method = 0,
        Use_Halo_Boost = False,
        zmin = 8,
        zmax = 1000,
        nz = 1500,
        nm = 500,
        sbh_range = [-3, 5],
        show_status = False,
        Use_Halo_Boost_Interp = True,
        Use_Halo_Boost_Extrap = True,
        Use_mdot_Interp = True,
        Use_LowM_Extrap_mdot = True,
        Use_HighM_Extrap_mdot = True,
        Break_OmB_Lost_Ratio = 1.1,
        Show_Extrap_MSG = False,
        Fix_mdot = False,
        Fixed_mdot = 10,
        Ignore_Evo = False):
    '''
    Find max fbh, don't think this is gonna be very helpful but let's see what happens
    '''
    Fbh_Range = np.array([1e-30, 0.99])
    LgFbh_Range = np.log10(Fbh_Range)
    
    def Radio_dif(LgF):
        r =  T_Radio(
            mc = mc,
            sbh = sbh,
            fbh0 = 10**LgF,
            Use_Ps_domain = False,
            OmB0 = OmB0,
            Use_EoR = Use_EoR,
            dmdm0_method = dmdm0_method,
            Use_Halo_Boost = Use_Halo_Boost,
            zmin = zmin,
            zmax = zmax,
            nz = nz,
            nm = nm,
            sbh_range = sbh_range,
            show_status = show_status,
            Use_Halo_Boost_Interp = Use_Halo_Boost_Interp,
            Use_Halo_Boost_Extrap = Use_Halo_Boost_Extrap,
            Use_mdot_Interp = Use_mdot_Interp,
            Use_LowM_Extrap_mdot = Use_LowM_Extrap_mdot,
            Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot,
            Break_OmB_Lost_Ratio = Break_OmB_Lost_Ratio,
            Show_Extrap_MSG = Show_Extrap_MSG,
            Fix_mdot = Fix_mdot,
            Fixed_mdot = Fixed_mdot,
            Ignore_Evo = Ignore_Evo)
        
        z = r['z']
        T = r['T']
        T0 = np.interp(x = 0, xp = z, fp = T)
        T_arcade = 0.467 # arcade excess at z0
        dif = T0-T_arcade
        
        return dif
    try:
        LgF_max = Solve(F = Radio_dif, Xmin = LgFbh_Range[0], Xmax = LgFbh_Range[1], Precision = 1e-3, show_status = show_status)
    except:
        dif1 = Radio_dif(LgFbh_Range[0])
        dif2 = Radio_dif(LgFbh_Range[1])
        warnings.warn('Solution not found, returning NaN. Info:')
        print('mc = ', mc, ', dif1 = ', dif1, ', dif2 = ', dif2)
        LgF_max = np.nan
    r = 10**LgF_max
    return r
