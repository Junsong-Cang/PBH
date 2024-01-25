'''
Main module for computing PBH accretion history,
written by JSC
'''

import numpy as np
import warnings, platform, ctypes
from Useful_Numbers import Cosmology as cosmo
from Useful_Numbers import Constants
from p21c_tools import *
try:
    from src.mdot_kernel import *
except:
    # for linked import
    from mdot_kernel import *

if platform.system() == 'Darwin':
    main_path = '/Users/cangtao/cloud/GitHub/PBH/'
else:
    # running on IHEP server
    main_path = '/home/dm/watson/work/Accreting_PBH/'

c_function_lib_file = main_path + 'src/c/mdot_interp.so'
c_function_lib = ctypes.CDLL(c_function_lib_file)

def Get_mdot_c(
        z = 10,
        m = 1000,
        OmB = cosmo['OmB'],
        Use_LowM_Extrap = True,
        Use_HighM_Extrap = True,
        Use_EoR = True,
        Use_Edd_Limit = True):
    '''
    c version of Get_mdot
    libs needs to be loaded externally for speed
    speed : 1 million calls per second, yes you heard that right
    '''
    
    if Use_LowM_Extrap:
        LowM_Extrap = 1
    else:
        LowM_Extrap = 0
    if Use_HighM_Extrap:
        HighM_Extrap = 1
    else:
        HighM_Extrap = 0
    if Use_EoR:
        EoR = 1
    else:
        EoR = 0
    if Use_Edd_Limit:
        Edd = 1
    else:
        Edd = 0

    # Load the shared library
    
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
    r = c_function(z, m, OmB, EoR, Edd, LowM_Extrap, HighM_Extrap)
    return r

def Get_Transfer_Data():
    TransferFile = main_path + 'data/Transfer_Function.h5'
    f = h5py.File(TransferFile, 'r')
    T_Array = f['T_Array'][:]
    Ek_axis_GeV = f['Axis/Kinetic_Energy_GeV'][:]
    Ek_axis = Ek_axis_GeV*1e9
    z_axis = f['Axis/z'][:]
    Hubble_axis = Hubble(z_axis)
    f.close()
    nz = len(z_axis)
    ne = len(Ek_axis)
    a = 1/(1+z_axis)
    lna = np.log(a)
    dlna_vec = (lna[0: nz-1] - lna[1:nz])
    dlna = np.sum(dlna_vec)/nz
    T_awap = T_Array/dlna
    T = np.empty((5, nz, ne, nz))

    for cid in np.arange(0,5):
        for zid_dep in np.arange(0, nz):
            for eid in np.arange(0, ne):
                for zid_inj in np.arange(0, nz):
                    T[cid, zid_dep, eid, zid_inj] = T_awap[zid_inj, eid, zid_dep, cid, 0]
    
    # tmp = T[3, 32, 10, 43]*dlna
    # print(tmp)

    TransferData = {
        'T' : T,
        'z_axis' : z_axis,
        'Ek_axis' : Ek_axis,
        'Hubble_axis' : Hubble_axis
        }
    return TransferData

def Get_Interp_Data():
    '''
    Get interpolation data for mdot, Boost Factor, Transfer function
    '''
    DataPath = main_path + 'data/'
    
    # mdot interp data, author : $main_path/mdot_main.py
    # contents:
    # mdot_vec : mdot without Edd Limit, index: [EoR, z, mbh]
    # z : z axis
    # m : m axis
    
    mdot_tab_file = DataPath + 'mdot_tab.npz'
    mdot_data = np.load(mdot_tab_file)

    # Boost Factor Table, author: Boost_Tab.py
    datafile = DataPath + 'Boost_Tab.npz'
    Boost_Factor_Tab = np.load(datafile)

    # Transfer data
    TransferData = Get_Transfer_Data()

    # Put them all in here
    r = {
        'mdot_data' : mdot_data,
        'Boost_Factor_Tab': Boost_Factor_Tab,
        'TransferData' : TransferData
        }
    return r

# Use as global variable for speed (avoid reload)
Interpolation_Table = Get_Interp_Data()

def Get_mdot(
        z = 10,
        m = 100,
        OmB = cosmo['OmB'],
        Use_C = True,
        Do_not_Interpolate = False,
        Use_LowM_Extrap = True,
        Use_HighM_Extrap = True,
        Use_EoR = True,
        Use_Edd_Limit = True,
        Show_Extrap_MSG = False):
    '''
    A fast interface for mdot,
    use interpolation for low masses to speed up
    Interpolation precision: to be tested
    Interpolation speed(calls per second 1 cpu, MBP M2 Ultra):
    Use_C = 0: 1850
    Use_C = 1: 1e6
    maximum difference between Use_C=0 and Use_C=1 is about 10^-7
    ---- inputs ----
    z : redshift
    m : bh mass in msun
    Do_not_Interpolate : abort interpolation, useful for testing
    Use_LowM_Extrap : use extrapolation for LowM overflow
    Use_HighM_Extrap : use extrapolation for HighM overflow
    '''
    Small = 1e-250
    OmB_LCDM = cosmo['OmB']
    # Set a floor for OmB_Ratio, sometimes in numerical integration OmB may transite from Small to -Small and this will lead to errors
    OmB_Ratio = OmB/OmB_LCDM

    if OmB_Ratio < Small:
        return 0
    
    if Use_C and not Do_not_Interpolate:
        r = Get_mdot_c(z = z, m = m, OmB = OmB, Use_LowM_Extrap = Use_LowM_Extrap, Use_HighM_Extrap = Use_HighM_Extrap, 
                       Use_EoR = Use_EoR, Use_Edd_Limit = Use_Edd_Limit)
        if r > -Small:
            return r
    
    zp = 1 + z
    
    # Allowed overflow fraction, will be useful for later validation
    Extrap_Redundancy = 1e6
    Axis_Redundancy = 1e-5
    
    # Read interp data
    zp_axis = 1 + Interpolation_Table['mdot_data']['z']
    m_axis = Interpolation_Table['mdot_data']['m']
    
    if zp > zp_axis[-1] and not Do_not_Interpolate:
        raise Exception('z not in interp range')
    
    if Use_EoR:
        mdot_tab = Interpolation_Table['mdot_data']['mdot_vec'][1,:,:]
    else:
        mdot_tab = Interpolation_Table['mdot_data']['mdot_vec'][0,:,:]

    if not Do_not_Interpolate and Within_Range(x = zp, x_array = zp_axis):

        mbh = m

        if mbh < m_axis[0] and Use_LowM_Extrap:
            OverFlow = m_axis[0]/mbh
            if OverFlow < Extrap_Redundancy:
                mbh = m_axis[0]*(1+Axis_Redundancy)
                # don't always send warning
                if Show_Extrap_MSG:
                    print('Use_LowM_Extrap triggered, input : [z, m] = [', '{0:.3f}'.format(z), ',',  '{0:.4e}'.format(m), ']',
                        ', interp MinM: ', '{0:.4e}'.format(m_axis[0]), ', Overflow :', '{0:.4e}'.format(OverFlow))
        
        if mbh > m_axis[-1] and Use_HighM_Extrap:
            OverFlow = mbh/m_axis[-1]
            if OverFlow < Extrap_Redundancy:
                mbh = m_axis[-1]*(1-Axis_Redundancy)
                if Show_Extrap_MSG:
                    print('Use_HighM_Extrap triggered, input : [z, m] = [', '{0:.3f}'.format(z), ',',  '{0:.4e}'.format(m), ']', 
                        ', interp MaxM: ', '{0:.4e}'.format(m_axis[-1]), ', Overflow :', '{0:.4e}'.format(OverFlow))

        if Within_Range(mbh, m_axis):
            r = Interp_2D(
                Tab = mdot_tab,
                x_axis = zp_axis,
                y_axis = m_axis,
                x_target = zp,
                y_target = mbh,
                Use_Log_X = True,
                Use_Log_Y = True,
                Use_Log_Z = True
            )

            if Use_Edd_Limit:
                r = min(OmB_Ratio * r, 10)
            else:
                r = OmB_Ratio * r
            return r

    if not Do_not_Interpolate:
        if Use_HighM_Extrap or Use_LowM_Extrap:
            warnings.warn('Extrap failed to catch overflow within allowed range, I am handling this to mdot_kernel().:')
        else:
            warnings.warn('No match found in interpolation table, diverting to mdot_kernel(). Possible cause: Extrapolation deactivated, severe overflow')

    # If not in range and not Extrap then following lines will be executed
    # the inputs bellow must conform with interpolation table
    r = mdot_kernel(
        z = z,
        M_PBH = m,
        Omega_b = OmB,
        Use_halo = 1,
        Use_feedback = 0,
        Use_EoR = Use_EoR
        )
    if Use_Edd_Limit:
        r = min(r, 10)
    return r

def Boost_Factor_Lite(z = 10,
                 OmB = cosmo['OmB'],
                 Use_EoR = True,
                 nmh = 1000,
                 show_status = 0,
                 ncpu = 1
                 ):
    '''
    Get the Boost Factor assuming mass growth rate scales with baryon density,
    should give BoostFactor(Use_Edd_Limit = 0) results
    '''
    
    Is_Scalar(x = z, ErrorMethod = 2)
    OmC = cosmo['OmC']
    h = cosmo['h']
    OmM = OmB + OmC
    RhoCr = 2.775e11 * h**2 # msun/Mpc^3
    RhoM_cmv = RhoCr*OmM
    
    zp = 1+z
    xe, Tk = LCDM_HyRec(z = z, Use_EoR = Use_EoR)
    Mmin = 1.3e3 * (10/zp)**1.5 * Tk**1.5
    t1 = TimeNow()
    m, dndm = HMF(
        z = z,
        model = 1, # ST model
        Mmin = Mmin,
        Mmax = 1e22, # should be enough
        nm = nmh,
        POWER_SPECTRUM = 0,
        Use_Interp = 1
    )
    if show_status:
        print('HMF time = ', TimeNow() - t1)

    # Collapse fraction
    fcoll = np.trapz(x = m, y = m*dndm) / RhoM_cmv

    # Now do the in halo part
    def Halo_Profile_Integrator(mh = 1e4):

        ProFile = HaloProfile(z = z, mh = mh, OmB = OmB, nr = 1000)
        r = ProFile[0,:]
        RhoM = ProFile[1,:]
        RhoC = ProFile[2,:]
        RhoB = ProFile[3,:]
        
        fun = RhoB * RhoC * r**2
        result = np.trapz(x = r, y = fun)
        
        # result has dimension M^2/L^3 in msun^2/pc^3, now convert to msun^2/Mpc^3
        
        result = result * 1e18
        
        return result

    # okey dokey, now do hmf
    if ncpu == 1:
        nm = len(m)
        kernel_0 = np.linspace(0, 1, nm)
        for mid in np.arange(0, nm):
            kernel_0[mid] = Halo_Profile_Integrator(mh = m[mid])
        kernel2 = kernel_0 * dndm
    else:
        t1 = TimeNow()
        kernel_0 = Parallel(n_jobs=ncpu)(delayed(Halo_Profile_Integrator)(mh = x) for x in m)
        kernel2 = kernel_0 * dndm
        if show_status:
            print('Total time = ', TimeNow() - t1)
    
    hmf_integration = np.trapz(x = m, y = kernel2) # dimension : (msun/Mpc^3)^2, or rho^2
    Halo_part = 4 * np.pi * hmf_integration/( OmB * OmC * RhoCr**2 * zp**3)
    IGM_part = (1-fcoll)**2
    
    Total_Boost = Halo_part + IGM_part

    # return individual parts for possible interpolation needs
    r = [Total_Boost, Halo_part, IGM_part]
    r = max(Total_Boost, 1)

    return r

def BoostFactor(
        z = 10,
        m = 0.1,
        OmB = cosmo['OmB'],
        Use_EoR = True,
        nmh = 1000,
        show_status = 0,
        ncpu = 1,
        Use_Edd_Limit = True,
        Use_Interp = True,
        Use_Extrap = True,
        Use_mdot_Interp = True,
        Use_LowM_Extrap_mdot = True,
        Use_HighM_Extrap_mdot = True
    ):
    '''
    Boost in mass growth rate caused by DM halos,
    Note that this is the globally averaged boost, growth rate for individual BH can scatter depending on the environment
    ----inputs----
    z : z
    m : BH mass in msun
    Use_EoR : Use EoR in 1. Minimum Virial halo mass; 2. mdot
    nmh : number of halo mass samples in HMF
    Use_Halo_Profile_Interp : Use interpolation, probably will be deappreciated bacause this does not make codes go faster
    show_status : show status
    ncpu : ncpu, must be set to 1 if using external MPI
    Use_Edd_Limit : Whether to use Eddington Limit, if not then BH mass growth rate scales with baryon density
    Use_Interp : Use interpolation
    Use_OmB_Extrap : Use extrapolation in OmB, happens in experimentally excluded regime so this should be safe for param scan
    ----outputs----
    Halo Boost Factor
    '''
    
    Is_Scalar(x = z, ErrorMethod = 2)
    zp = 1+z
    
    if Use_Interp:
        Boost_Table = Interpolation_Table['Boost_Factor_Tab']
        Boost_Tab = Boost_Table['Boost_Tab']
        OmB_axis = Boost_Table['OmB_vec']
        m_axis = Boost_Table['m_vec']
        zp_axis = Boost_Table['zp_vec']
        # check ranges
        nz, nm, nb = len(zp_axis), len(m_axis), len(OmB_axis)
        m_ok = Within_Range(m, m_axis)
        z_ok = Within_Range(zp, zp_axis)
        OmB_ok = Within_Range(OmB, OmB_axis)
        
        Tab = np.empty((nb, nz, nm))
        if Use_EoR:
            Tab[:,:,:] = Boost_Tab[1,:,:,:]
        else:
            Tab[:,:,:] = Boost_Tab[0,:,:,:]
        # Let's do it in log axis, the minimum of Tab is 1 so we can do Tab in log too
        lr = Interp_3D(
            Tab = np.log10(Tab),
            x_axis = np.log10(OmB_axis),
            y_axis = zp_axis - 1,
            z_axis = np.log10(m_axis),
            x = np.log10(OmB),
            y = z,
            z = np.log10(m),
            Use_Extrap = True
        )
        result = 10**lr
        
        # decide whether to use this result

        if m_ok and z_ok and OmB_ok:
            return result
            
        # If the code execute following block it means overflow detected, but there are physically motivated ways to get around it        
        # Only z<Zmin will get diverted to slow mode, more explanations will be added
        if zp < zp_axis[0]:
            warnings.warn('z is too small for interpolation, the rest would be sloooooooow')
        else:
            if not Use_Extrap:
                warnings.warn('Interpolation overflow detected, consider enabling Use_Extrap to use extrapolation (safe for params scan)')
            else:
                return result
        
    OmC = cosmo['OmC']
    h = cosmo['h']
    OmM = OmB + OmC
    RhoCr = 2.775e11 * h**2 # msun/Mpc^3
    RhoM_cmv = RhoCr*OmM
    RhoC_cmv = RhoCr*OmC
    RhoB_avg = RhoCr*OmB*(1+z)**3
    
    mdot_ideal = Get_mdot(
        z = z, 
        m = m, 
        OmB = OmB, 
        Use_EoR = Use_EoR, 
        Do_not_Interpolate = not Use_mdot_Interp,
        Use_LowM_Extrap = Use_LowM_Extrap_mdot,
        Use_HighM_Extrap = Use_HighM_Extrap_mdot,
        Use_Edd_Limit = False)
    
    def Find_mdot_effective(RhoB):
        # RhoB : msun/Mpc^3
        ratio = RhoB/RhoB_avg
        mdot_scaled = ratio * mdot_ideal
        if not Use_Edd_Limit:
            return mdot_scaled

        subtract_factor = (mdot_scaled-10) * np.heaviside(mdot_scaled-10, 0)
        mdot_edd_limited = mdot_scaled - subtract_factor # this should ensure a max of 10
        # print('mdot_max = ', np.max(mdot_edd_limited))
        # now substract radiation
        # f_esc = min(0.1, m_dot) # escape fraction
        f_esc_bright = 0.1 * np.heaviside(mdot_edd_limited- 0.1, 0)
        f_esc_dim = mdot_edd_limited * np.heaviside(0.1 - mdot_edd_limited, 0)
        f_esc = f_esc_bright + f_esc_dim
        mdot_effective = mdot_edd_limited*(1-f_esc)
        # print('mdot_effective = ', np.max(mdot_effective))
        
        return mdot_effective
        
    xe, Tk = LCDM_HyRec(z = z, Use_EoR = Use_EoR)
    Mmin = 1.3e3 * (10/zp)**1.5 * Tk**1.5
    m, dndm = HMF(
        z = z,
        model = 1, # ST model
        Mmin = Mmin,
        Mmax = 1e22, # should be enough
        nm = nmh,
        POWER_SPECTRUM = 0,
        Use_Interp = 1
    )
    
    # ----HMG----
    mdot_HMG = Find_mdot_effective(RhoB_avg)
    B_HMG = mdot_HMG * RhoC_cmv
    
    # ----IGM----
    fcoll = np.trapz(x = m, y = m*dndm) / RhoM_cmv
    RhoB_IGM = RhoB_avg * (1 - fcoll)
    mdot_IGM = Find_mdot_effective(RhoB_IGM)
    B_IGM = mdot_IGM * (1 - fcoll)*RhoC_cmv

    # ----Halo----
    def Halo_I_Factor(mh):
        Profile = HaloProfile(
            z = z,
            mh = mh,
            OmB = OmB,
            nr = 1000,
            map_nx = 100,
            map_precision = 1e-4,
            mass_error = 0.05
            )
        # All in Mpc unit
        r = Profile[0,:] / 1e6
        RhoM = Profile[1,:] * 1e18
        RhoC = Profile[2,:] * 1e18
        RhoB = Profile[3,:] * 1e18
        mdot_Halo = Find_mdot_effective(RhoB)
        fun = r**2 * mdot_Halo * RhoC
        result = np.trapz(x = r, y = fun)
        return result
    
    # okey dokey, now do hmf
    if ncpu == 1:
        IH = np.linspace(0, 1, nmh)
        for mid in np.arange(0, nmh):
            IH[mid] = Halo_I_Factor(mh = m[mid])
        kernel =  4 * np.pi * dndm * IH
    else:
        t1 = TimeNow()
        IH = Parallel(n_jobs=ncpu)(delayed(Halo_I_Factor)(mh = x) for x in m)
        kernel =  4 * np.pi * dndm * IH
    
    B_Halo = np.trapz(x = m, y = kernel)
    
    result = (B_IGM + B_Halo)/B_HMG
    result = max(result, 1)
    
    return result

def dmdx(
        x = -2,
        m = 100,
        OmB = cosmo['OmB'],
        Use_EoR = True,
        Fix_mdot = False,
        Fixed_mdot = 10,
        Use_Halo_Boost = False,
        Use_Halo_Boost_Interp = True,
        Use_Halo_Boost_Extrap = True,
        Use_mdot_Interp = True,
        Use_LowM_Extrap_mdot = True,
        Use_HighM_Extrap_mdot = True,
        Show_Extrap_MSG = False
        ):
    '''
    Get net bh mass growth rate dm/dlna and luminosity in J/s
    ----inputs----
    x : lna
    m : mbh in msun
    Use_EoR : Use EoR
    Fix_mdot : debug option, whether to fix mdot
    Fixed_mdot : debug option, fixed value for mdot
    Use_Halo_Boost : Use_Halo_Boost, needs to be handled carefully as it indicates a environment-dependent scatter in growth rate
    Use_Halo_Boost_Interp : Use interpolation in BoostFactor
    ----outputs----
    dmdx
    Lx : x-ray bolometric luminosity in J/s
    LR : Radio luminosity at 1 GHz, in J/s/Hz
    '''

    msun = Constants['msun']
    c = Constants['c']
    aR = 0.62 # Radio spectra index, L_R \propto v^{-aR}

    a = np.exp(x)
    z = 1/a - 1
    mdot_edd = 1.44e14 * m / msun # edd rate in msun/s
    Ledd = mdot_edd * msun * c**2 # J/s

    if Fix_mdot:
        m_dot = Fixed_mdot
    else:
        m_dot = Get_mdot(
            m = m, 
            z = z, 
            OmB = OmB,
            Use_EoR = Use_EoR,
            Do_not_Interpolate = not Use_mdot_Interp,
            Use_LowM_Extrap = Use_LowM_Extrap_mdot,
            Use_HighM_Extrap = Use_HighM_Extrap_mdot,
            Show_Extrap_MSG = Show_Extrap_MSG
            )

    f_esc = min(0.1, m_dot) # escape fraction
    m_dot_effective = m_dot * (1-f_esc)
    
    dmdt = m_dot_effective * mdot_edd # msun/s

    if Use_Halo_Boost:
        dmdt = dmdt * BoostFactor(
            z = z,
            m = m,
            OmB = OmB,
            Use_EoR = Use_EoR,
            Use_Interp = Use_Halo_Boost_Interp,
            Use_Extrap = Use_Halo_Boost_Extrap,
            Use_mdot_Interp = Use_mdot_Interp,
            Use_LowM_Extrap_mdot = Use_LowM_Extrap_mdot,
            Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot
            )
    
    if m_dot <= 0.02:
        Lbh_ratio = 2.3 * np.log10(m_dot) + 1.1
    else:
        Lbh_ratio = 0.25 * np.log10(m_dot) - 2.4
    
    Lx = 10**Lbh_ratio * Ledd

    # Now get radio,ref: zzh2023
    Lx_erg = Lx*1e7
    LgL5 = 0.6 * np.log10(Lx_erg) + 0.78 * np.log10(m) + 7.33
    Lb5 = 10.0**LgL5 / 1.0e7 # f*LR @ 5GHz, in J/s, f is frequency
    f5 = 5.0e9
    L5 = Lb5/f5 # dEdfdt @ 5GHz, in J/s/Hz
    LR = L5*5.0**aR

    H = Hubble(z)
    r = dmdt/H, Lx, LR
        
    return r

def dOmBdx(
        dm_dx = 0,
        fbh0 = 1e-4,
        m0 = 1e6
        ):
    '''
    Get dOmB/dlna
    '''
    OmC = cosmo['OmC']
    r = - dm_dx * fbh0 * OmC / m0
    return r

def z2x(z):
    '''
    Convert z to x=lna
    '''
    a = 1/(1+z)
    r = np.log(a)
    return r

def x2z(x):
    '''
    Convert lna to z
    '''
    a = np.exp(x)
    r = 1/a - 1
    return r

def Get_Evo_Lite(
        fbh0 = 1e-5,
        m0 = 1e3,
        OmB0 = cosmo['OmB'],
        zmin = 8,
        zmax = 1000,
        nz = 2000,
        Use_Halo_Boost = False,
        Use_EoR = True,
        Use_OmB_FeedBack = True,
        Fix_mdot = False,
        Fixed_mdot = 10,
        Use_Halo_Boost_Interp = True,
        Use_Halo_Boost_Extrap = True,
        Use_mdot_Interp = True,
        Use_LowM_Extrap_mdot = True,
        Use_HighM_Extrap_mdot = True,
        Show_Extrap_MSG = False,
        Ignore_Evo = False
    ):
    '''
    Find mass and OmB evolution for delta mass function
    ----inputs----
    fbh0 : initial fbh
    m0 : initial mass at zmax
    OmB0 : Initial OmegaB
    zmin : minimum of redshift integration range
    zmax : maximax of redshift integration range
    nz : z timesteps
    Use_OmB_FeedBack : use OmegaB feedback
    Fix_mdot : test feature, whether or not to fix mdot
    Fixed_mdot : fixed mdot
    Use_Halo_Boost : Use_Halo_Boost
    Ignore_Evo : do not evolve mass, use Lbh only
    ----outputs----
    various ratios and luminosity
    '''
    
    Small = 1e-280
    x1 = z2x(zmax)
    x2 = z2x(zmin)
    dx = (x2 - x1)/(nz - 1)
    x_vec = np.linspace(x1, x2, nz)
    z_vec = x2z(x_vec)
    
    m_vec = np.linspace(0, 1, nz)
    OmB_vec = np.linspace(0, 1, nz)
    dmdx_vec = np.linspace(0, 1, nz)
    Lbh_vec = np.linspace(0, 1, nz)
    LR_vec = np.linspace(0, 1, nz)

    m_vec[0] = m0
    OmB_vec[0] = OmB0
    dmdx_vec[0] = 0
    Lbh_vec[0] = 0
    LR_vec[0] = 0
    
    # Start evolving
    for xid in np.arange(1, nz):
        xid_prev = xid-1

        x_prev = x_vec[xid_prev]
        OmB_prev = OmB_vec[xid_prev]
        m_prev = m_vec[xid_prev]
        if OmB_prev < Small:
            dm_dx, Lbh, LR = 0.0, 0.0, 0.0
        else:
            dm_dx, Lbh, LR = dmdx(
                x = x_prev,
                m = m_prev,
                OmB = OmB_prev,
                Use_EoR = Use_EoR,
                Fix_mdot = Fix_mdot,
                Fixed_mdot = Fixed_mdot,
                Use_Halo_Boost = Use_Halo_Boost,
                Use_Halo_Boost_Interp = Use_Halo_Boost_Interp,
                Use_Halo_Boost_Extrap = Use_Halo_Boost_Extrap,
                Use_mdot_Interp = Use_mdot_Interp,
                Use_LowM_Extrap_mdot = Use_LowM_Extrap_mdot,
                Show_Extrap_MSG = Show_Extrap_MSG,
                Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot)
        if Ignore_Evo:
            dm_dx = 0
        dmdx_vec[xid] = dm_dx
        Lbh_vec[xid] = Lbh
        LR_vec[xid] = LR

        if Use_OmB_FeedBack:
            dOmB_dx = dOmBdx(
                dm_dx = dm_dx,
                fbh0 = fbh0,
                m0 = m0
            )
        else:
            # this is a test feature anyway so don't worry about fix_mdot
            dOmB_dx = 0

        dm = dm_dx * dx
        dOmB = dOmB_dx * dx
        
        m_vec[xid] = m_prev + dm
        OmB_vec[xid] = max(OmB_prev + dOmB, Small / 2)
    
    # Get dEdVdt
    OmC = cosmo['OmC']
    h = cosmo['h']
    RhoCr = 2.775e11 * h**2
    RhoC = RhoCr * OmC *(1+z_vec)**3 # msun/Mpc^3
    nbh = fbh0 * RhoC/m0 # Mpc^3
    
    Q = Constants['Q']
    Mpc = Constants['Mpc']
    Mpc_cm = Mpc*100
    
    Lbh_vec_eV = Lbh_vec/Q
    dEdVdt = nbh * Lbh_vec
    # Convert from J/s/Mpc^3 to eV/s/cm^3:
    dEdVdt = dEdVdt / Q / Mpc_cm**3
    
    # Now get radio
    Radio_EMS = nbh*LR_vec/(1+z_vec)**3 / Mpc**3
    
    r = {
        'z' : z_vec, 
        'm_ratio' : m_vec/m0, 
        'OmB_ratio' : OmB_vec/OmB0, 
        'fbh_ratio' : m_vec/m0, 
        'dmdx_vec' : dmdx_vec, 
        'x_vec' : x_vec, 
        'dEdVdt' : dEdVdt, # x-ray
        'Lbh_vec_eV' : Lbh_vec_eV, # x-ray
        'Radio_EMS' : Radio_EMS} # COMOVING Radio emissivity @ 1GHz, unit:J/s/m^3/Hz

    return r

def Find_fbhmax_lite(
        m0 = 1e3,
        OmB0 = cosmo['OmB'],
        OmB_Lost_Ratio = 0.3,
        Use_EoR = True,
        Use_Halo_Boost = False,
        zmin = 8,
        zmax = 1000,
        nz = 5000,
        lf_min = -30,
        Precision = 1e-2,
        Use_Halo_Boost_Interp = True,
        Use_Halo_Boost_Extrap = True,
        Use_mdot_Interp = True,
        Use_LowM_Extrap_mdot = True,
        Use_HighM_Extrap_mdot = True,
        Show_Extrap_MSG = False
        ):
    '''
    Find fbh limit for monochromatic BHs
    ----inputs----
    m0 : initial mass in msun
    OmB0 : Initial OmB
    OmB_Lost_Ratio : 
        maximumly allowed OmB lost ratio, default value of 0.3 is motivated by DES Low-Z results
        DOI: 10.1103/PhysRevD.105.023520
    Use_EoR : Use EoR
    Use_Halo_Boost : Use Halo Boost
    zmin : minimum evolve redshift
    zmax : maximum evolve redshift
    nz : number of z steps to use, log-distributed between [zmin, zmax]
    lf_min : fbh search minimum in log
    Precision : solution precision
    '''
    
    def omb_dif(LgFbh0):
        '''
        Find whether fraction fbh0 under-accret omb
        '''
        r = Get_Evo_Lite(
            fbh0 = 10**LgFbh0,
            m0 = m0,
            OmB0 = OmB0,
            zmin = zmin,
            zmax = zmax,
            nz = nz,
            Use_OmB_FeedBack = True,
            Use_EoR = Use_EoR,
            Use_Halo_Boost = Use_Halo_Boost,
            Use_Halo_Boost_Interp = Use_Halo_Boost_Interp,
            Use_Halo_Boost_Extrap = Use_Halo_Boost_Extrap,
            Use_mdot_Interp = Use_mdot_Interp,
            Show_Extrap_MSG = Show_Extrap_MSG,
            Use_LowM_Extrap_mdot = Use_LowM_Extrap_mdot,
            Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot
        )
        ob_acc_ratio = 1 - np.min(r['OmB_ratio'])
        r = ob_acc_ratio - OmB_Lost_Ratio
        return r
    
    # For small m0 (accretion inefficient) Solve might fail
    
    try:
        lg_fbh0_max = Solve(
            F = omb_dif,
            Xmin = lf_min,
            Xmax = 0,
            Precision = Precision)
    except:
        lg_fbh0_max = 0
    
    r = 10**lg_fbh0_max

    return r

def Log_Normal(m = 1e2, mc = 1e1, sbh = 0.5):
    '''
    Normalised log-normal profile dfbh/dlnm
    '''
    lnm = np.log(m)
    lnmc = np.log(mc)
    chi2 = (lnm - lnmc)**2/(2*sbh**2)
    # prefix = 1/(np.sqrt(2 * np.pi)*sbh*m)
    prefix = 1/(np.sqrt(2 * np.pi)*sbh)
    r = prefix*np.exp(-chi2)
    return r

def Get_Evo(
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
        Use_HZ = False,
        nz = 1000,
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
    '''
    Get evolution of PBH mass function and OmB
    Phi is defined as drho/dlnm/rho0

    ----inputs----

    mc : center mass in LN MF
    sbh : sigma_bh in LN MF
    fbh0 : initial fbh

    Use_Ps_domain : compute evolution using PBH profile generated by Ps params [PsA, kc, SigmaPs],
                    If True then [mc, sbh, fbh0] will be ignored
                    Default setting gives [fbh, mc, SigmaBH] = [2e-2, 576, 0.294]
    PsA : Amplitude A for curvature power spectrum Ps
    kc : central k for Ps, in Mpc^-1
    SigmaPs : Sigma for Ps
    DeltaC : PBH collapse threshold

    OmB0 : initial OmegaB
    Use_EoR : Use EoR
    Use_Halo_Boost : Use Halo Boost Factor
    zmin : minimum z
    zmax : maximum z
    Use_HZ : Use Lbh at z > 1000, useful for recombination calculations
    nz : z time steps
    nm : mbh time steps
    sbh_range : how many sigmas to use in mbh range
    show_status : show status of evolution

    Use_Halo_Boost_Interp : Use Interpolation in Halo Boost
    Use_Halo_Boost_Extrap : Use Extrapolation in Halo Boost
    Use_mdot_Interp : Use Interpolation in mdot calculation
    Use_LowM_Extrap_mdot : Use LowM Extrapolation in mdot calculation
    Use_HighM_Extrap_mdot : Use HighM Extrapolation in mdot calculation
    Break_OmB_Lost_Ratio : set dmdt and Lbh to 0 when baryon lost ratio exceeds this, can speed up calculation for param scan, set to >1 to use no cut
    Fix_mdot : whether to fix mdot, useful for testing
    Fixed_mdot : fixed value for mdot
    Ignore_Evo : do not evolve mass, use Lbh only
    '''
    
    OmC = cosmo['OmC']
    h = cosmo['h']
    Q = Constants['Q']
    Mpc = Constants['Mpc']
    Mpc_cm = Mpc * 100.0

    RhoCr = 2.775e11 * h**2 # msun/Mpc^3
    RhoC_cmv = RhoCr*OmC
    
    Small = 1e-280
        
    # Get x array for time-evolution
    x1 = z2x(zmax)
    x2 = z2x(zmin)
    dx = (x2 - x1)/(nz - 1)
    x_vec = np.linspace(x1, x2, nz)
    z_vec = x2z(x_vec)
    
    # Get Phi0, either from log-normal profile or actual profile generated by Ps
    if not Use_Ps_domain:
        # Initialise Mass array
        yc = np.log(mc)
        y1 = yc + sbh*sbh_range[0]
        y2 = yc + sbh*sbh_range[1]
        lnm0 = np.linspace(y1, y2, nm)    
        m0_vec = np.exp(lnm0)
        Phi0 = Log_Normal(m = m0_vec, mc = mc, sbh = sbh)
        f0 = fbh0
    else:
        
        # Get fbh0 and profile shape
        PBH_Profile_Raw = Phi_Profile(
            A = PsA,
            Sigma = SigmaPs,
            kc = kc,
            DeltaC = DeltaC
        )
        
        f0 = PBH_Profile_Raw['fbh']
        SigmaBH = PBH_Profile_Raw['Width']

        if (f0 > 1e-60) and (not np.isnan(SigmaBH)):
            
            Mc = PBH_Profile_Raw['mc']
            # Note that Raw profile is normalised to fbh!
            Phi0_axis = PBH_Profile_Raw['Phi_vec']/f0
            m_axis = PBH_Profile_Raw['m_vec']

            # Possible alternative : use analytic profile from [mc, sbh] params
            yc = np.log(Mc)
            y1 = yc + SigmaBH * sbh_range[0]
            y2 = yc + SigmaBH * sbh_range[1]
            lnm0 = np.linspace(y1, y2, nm)
            m0_vec = np.exp(lnm0)
            Phi0 = np.linspace(0,1,nm)

            for idx in np.arange(0,nm):
                m_ = m0_vec[idx]
                if Within_Range(x = m_, x_array = m_axis):
                    Phi0[idx] = np.interp(
                        x = np.log10(m_),
                        xp = np.log10(m_axis),
                        fp = Phi0_axis
                        )
                else:
                    # nanana, don't forget to normalise fbh
                    Phi0[idx] = (1/f0)*Phi_Kernel(
                        A = PsA,
                        Sigma = SigmaPs,
                        kc = kc,
                        m = m_,
                        DeltaC = DeltaC
                    )
        else:
            warnings.warn('Input Ps generates no PBH (fbh<1e-60), I am returning 0')
            # Get some Fake Profile
            f0 = 0
            SigmaBH = 1
            Mc = k2m(kc)
            LgMc = np.log10(Mc)
            
            Phi0_axis = np.zeros(nm)
            m_axis = np.logspace(LgMc - 2, LgMc + 2, nm)
            yc = np.log(Mc)
            y1 = yc + SigmaBH * sbh_range[0]
            y2 = yc + SigmaBH * sbh_range[1]
            lnm0 = np.linspace(y1, y2, nm)
            m0_vec = np.exp(lnm0)
            Phi0 = np.zeros(nm)
            
        # Ok Phi0 now obtained
    
    # Create some arrays
    omb_vec = np.linspace(0, 1, nz)
    m_vec = np.empty((nz, nm))
    m_ratio = np.empty((nz, nm))
    Lbh_vec = np.zeros((nz, nm))
    Phi_vec = np.empty((nz, nm))
    Phi_here = np.linspace(0, 1, nm)
    m_now = np.linspace(0, 1, nm)
    dEdVdt_mono = np.linspace(0, 1, nm)
    EMS_mono = np.linspace(0, 1, nm)
    fbh_ratio = np.linspace(0, 1, nz)
    dEdVdt = np.linspace(0, 1, nz)
    Radio_EMS = np.zeros(nz)
    
    # Initial Condition
    m_vec[0,:] = m0_vec[:]
    omb_vec[0] = OmB0
    fbh_ratio[0] = 1

    if show_status:
        t1 = TimeNow()
        print(' Mass range: [', m0_vec[0], ', ', m0_vec[-1])

    RhoBH0_cmv = f0 * RhoC_cmv
    Phi_vec[0,:] = Phi0[:]
    
    for zid in np.arange(1, nz):
        if show_status:
            print('Get_Evo status : ', zid/nz)
        zid_prev = zid-1
        x_prev = x_vec[zid_prev]
        omb_prev = omb_vec[zid_prev]
        OmB_Lost_Fraction = (OmB0 - omb_prev)/OmB0
        
        # Find mass evo
        for mid in np.arange(0, nm):
            try:
                if (omb_prev < Small) or (f0 < Small):
                    # Numerical cut, this means there is no baryon or BH left
                    dm_dx, Lbh, LR = 0.0, 0.0, 0.0
                elif OmB_Lost_Fraction > 1.05*Break_OmB_Lost_Ratio:
                    # Economical cut
                    dm_dx, Lbh, LR = 0.0, 0.0, 0.0
                else:
                    dm_dx, Lbh, LR = dmdx(
                        x = x_prev, 
                        m = m_vec[zid_prev, mid], 
                        OmB = omb_prev, 
                        Use_EoR = Use_EoR, 
                        Use_Halo_Boost = Use_Halo_Boost,
                        Use_Halo_Boost_Interp = Use_Halo_Boost_Interp,
                        Use_Halo_Boost_Extrap = Use_Halo_Boost_Extrap,
                        Use_mdot_Interp = Use_mdot_Interp,
                        Use_LowM_Extrap_mdot = Use_LowM_Extrap_mdot,
                        Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot,
                        Fix_mdot = Fix_mdot,
                        Fixed_mdot = Fixed_mdot,
                        Show_Extrap_MSG = Show_Extrap_MSG
                        )
                    if Ignore_Evo:
                        dm_dx = 0
            except:
                print('----crash iminent----')
                print('z = ', x2z(x_prev), ', m = ', m_vec[zid_prev, mid], ', OmB_ratio = ', omb_prev/OmB0)
                raise Exception('dmdx crashed')
            
            dm = dm_dx*dx
            m_vec[zid, mid] = m_vec[zid_prev, mid] + dm
            m_ratio[zid, mid] = m_vec[zid, mid]/m0_vec[mid]
            Lbh_vec[zid, mid] = Lbh
            
            # Find dEdVdt for mono
            nbh = RhoBH0_cmv * ((1+z_vec[zid])**3)/m0_vec[mid] # Mpc^-3
            swap = nbh * Lbh # J/s/Mpc^3
            swap = swap/Q/Mpc_cm**3 # eV/s/cm^3
            dEdVdt_mono[mid] = swap

            # Find radio EMS for mono
            nbh_cmv = nbh/(1+z_vec[zid])**3/Mpc**3 # comoving bh number density, m^-3
            EMS_mono[mid] = nbh_cmv * LR # comoving radio EMS, J/s/Hz/m^3

        dEdVdt[zid] = np.trapz(x = lnm0, y = dEdVdt_mono * Phi0)
        Radio_EMS[zid] = np.trapz(x = lnm0, y = EMS_mono * Phi0)

        # Find MF evo
        m_now[:] = m_vec[zid,:]
        lnm = np.log(m_now)
        dlnm0_dlnm = Get_dydx(x = lnm, y = lnm0, Use_log = False, method = dmdm0_method)
        if np.min(dlnm0_dlnm) < -Small:
            raise Exception('Mass overlap detected, smaller BH can grow faster?')
        Phi_here = Phi0 * dlnm0_dlnm * m_now/m0_vec # new profile
        if np.min(Phi_here) < -Small:
            raise Exception('Phi is negative, this is not supposed to happen.')
        if True in np.isnan(Phi_here):
            warnings.warn('Nan found in Phi_here')
            print('Phi_here = ', Phi_here)
            print('Phi0 = ', Phi0)
            print('dlnm0_dlnm = ', dlnm0_dlnm)
            print('m_now = ', m_now)
            print('m0_vec = ', m0_vec)
            print('omb_prev = ', omb_prev)
            
            raise Exception('Nan detected')

        Phi_vec[zid,:] = Phi_here[:]

        # Find OmB
        fbh_ratio[zid] = np.trapz(x = lnm, y = Phi_here) # fbh/fbh0
        fbh_prev = f0 * fbh_ratio[zid_prev]
        fbh_now = f0 * fbh_ratio[zid]
        # conservation of matter dictates that 
        # OmB2 = OmB1 + (fbh1 - fbh2)*OmC
        omb_now = omb_prev + (fbh_prev - fbh_now)*OmC
        omb_vec[zid] = max(omb_now, Small/2)
    
    dEdVdt[0] = dEdVdt[1]
    Lbh_vec[0,:] = Lbh_vec[1,:]
    z_axis_max = np.max(Interpolation_Table['mdot_data']['z'])
    if Use_HZ and (zmax < z_axis_max):
        
        nzh = 60
        
        z_vec_hz = np.linspace(z_axis_max - 0.001, zmax + 0.01, nzh)
        omb_vec_hz = OmB0*np.ones(nzh)
        fbh_ratio_hz = np.ones(nzh)
        m_ratio_hz = np.ones((nzh, nm))
    
        Phi_vec_hz = np.zeros((nzh, nm))
        m_vec_hz = np.zeros((nzh, nm))
        
        for zid in np.arange(0, nzh):
            Phi_vec_hz[zid,:] = Phi0[:]
            m_vec_hz[zid,:] = m0_vec[:]
        
        # All I really want is Lbh and perhaps dEdVdt
        dEdVdt_hz = np.zeros(nzh)
        Radio_EMS_hz = np.zeros(nzh)
        
        Lbh_vec_hz = np.zeros((nzh, nm))

        for zid in np.arange(0, nzh):

            z_ = z_vec_hz[zid]
            x_ = z2x(z_)
                
            for mid in np.arange(0, nm):
                
                m_ = m0_vec[mid]
                
                dm_dx, Lbh, LR = dmdx(
                        x = x_, 
                        m = m_,
                        OmB = OmB0,
                        Use_EoR = Use_EoR, 
                        Use_Halo_Boost = Use_Halo_Boost,
                        Use_Halo_Boost_Interp = Use_Halo_Boost_Interp,
                        Use_Halo_Boost_Extrap = Use_Halo_Boost_Extrap,
                        Use_mdot_Interp = Use_mdot_Interp,
                        Use_LowM_Extrap_mdot = Use_LowM_Extrap_mdot,
                        Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot,
                        Fix_mdot = Fix_mdot,
                        Fixed_mdot = Fixed_mdot,
                        Show_Extrap_MSG = Show_Extrap_MSG
                        )
                
                Lbh_vec_hz[zid, mid] = Lbh
                nbh = RhoBH0_cmv * ((1+z_vec_hz[zid])**3)/m0_vec[mid] # Mpc^-3
                swap = nbh * Lbh # J/s/Mpc^3
                swap = swap/Q/Mpc_cm**3 # eV/s/cm^3
                dEdVdt_mono[mid] = swap
                nbh_cmv = nbh/(1+z_vec[zid])**3/Mpc**3

                EMS_mono[mid] = nbh_cmv * LR

            dEdVdt_hz[zid] = np.trapz(x = lnm0, y = dEdVdt_mono * Phi0)
            Radio_EMS_hz[zid] = np.trapz(x = lnm0, y = EMS_mono * Phi0)
        
        nz_full = nz+nzh

        z_vec_full = np.linspace(0, 10, nz_full)
        omb_vec_full = OmB0*np.ones(nz_full)
        fbh_ratio_full = np.ones(nz_full)
        m_ratio_full = np.ones((nz_full, nm))
    
        Phi_vec_full = np.zeros((nz_full, nm))
        m_vec_full = np.zeros((nz_full, nm))
        dEdVdt_full = np.zeros(nz_full)
        Radio_EMS_full = np.zeros(nz_full)
        Lbh_vec_full = np.zeros((nz_full, nm))

        # copying HZ
        z_vec_full[0:nzh] = z_vec_hz[0:nzh]
        omb_vec_full[0:nzh] = omb_vec_hz[0:nzh]
        fbh_ratio_full[0:nzh] = fbh_ratio_hz[0:nzh]
        dEdVdt_full[0:nzh] = dEdVdt_hz[0:nzh]
        Radio_EMS_full[0:nzh] = Radio_EMS_hz[0:nzh]
        m_ratio_full[0:nzh,:] = m_ratio_hz[0:nzh,:]
        Phi_vec_full[0:nzh,:] = Phi_vec_hz[0:nzh,:]
        m_vec_full[0:nzh,:] = m_vec_hz[0:nzh,:]
        Lbh_vec_full[0:nzh,:] = Lbh_vec_hz[0:nzh,:]

        # copying evolved table
        z_vec_full[nzh:nz_full] = z_vec[0:nz]
        omb_vec_full[nzh:nz_full] = omb_vec[0:nz]
        fbh_ratio_full[nzh:nz_full] = fbh_ratio[0:nz]
        dEdVdt_full[nzh:nz_full] = dEdVdt[0:nz]
        Radio_EMS_full[nzh:nz_full] = Radio_EMS[0:nz]
        m_ratio_full[nzh:nz_full,:] = m_ratio[0:nz,:]
        Phi_vec_full[nzh:nz_full,:] = Phi_vec[0:nz,:]
        m_vec_full[nzh:nz_full,:] = m_vec[0:nz,:]
        Lbh_vec_full[nzh:nz_full,:] = Lbh_vec[0:nz,:]

        # now rename them
        z_vec = z_vec_full
        omb_vec = omb_vec_full
        fbh_ratio = fbh_ratio_full
        dEdVdt = dEdVdt_full
        Radio_EMS = Radio_EMS_full
        m_ratio = m_ratio_full
        Phi_vec = Phi_vec_full
        m_vec = m_vec_full
        Lbh_vec = Lbh_vec_full
    
    Lbh_vec_eV = Lbh_vec/Q
    
    if show_status:
        t_used = TimeNow() - t1
        print('Get_Evo complete, time used: ', t_used)

    r = {
        'MF' : Phi_vec,
        'm_vec' : m_vec,
        'z_vec' : z_vec,
        'omb_ratio' : omb_vec/OmB0,
        'm0_vec' : m0_vec,
        'fbh_ratio' : fbh_ratio,
        'fbh_vec' : f0 * fbh_ratio,
        'm_ratio' : m_ratio,
        'dEdVdt' : dEdVdt,
        'Lbh_vec_eV' : Lbh_vec_eV, # x-ray bolometric luminosity, eV/s
        'Radio_EMS' : Radio_EMS # COMOVING Radio emissivity @ 1GHz, unit:J/s/m^3/Hz
        }

    return r

def Find_fbhmax(
        mc = 1e2,
        sbh = 0.5,
        OmB_Lost_Ratio = 0.3,
        Use_EoR = True,
        Use_Halo_Boost = False,
        sbh_range = [-3, 5],
        zmin = 8,
        zmax = 1000,
        nz = 1000,
        nm = 500,
        OmB0 = cosmo['OmB'],
        lf_min = -30,
        Precision = 1e-2,
        show_status = False,
        Use_Halo_Boost_Interp = True,
        Use_Halo_Boost_Extrap = True,
        Use_mdot_Interp = True,
        Use_LowM_Extrap_mdot = True,
        Use_HighM_Extrap_mdot = True,
        Show_Extrap_MSG = False,
        Fix_mdot = False,
        Fixed_mdot = 10
    ):
    '''
    Find maximum fbh0, takes about 14 Get_Evo calls
    ----inputs----
    OmB_Lost_Ratio : 
        maximumly allowed OmB lost ratio, default value of 0.3 is motivated by DES Low-Z results
        DOI: 10.1103/PhysRevD.105.023520
    Speed : 71s for nz = 1500, nm = 500
    '''
    
    def omb_dif(LgFbh0):
        r = Get_Evo(
            mc = mc,
            sbh = sbh,
            fbh0 = 10**LgFbh0,
            Use_Ps_domain = False,
            OmB0 = OmB0,
            Use_EoR = Use_EoR,
            Use_Halo_Boost = Use_Halo_Boost,
            zmin = zmin,
            zmax = zmax,
            nz = nz,
            nm = nm,
            sbh_range = sbh_range,
            show_status = False,
            Show_Extrap_MSG = Show_Extrap_MSG,
            Use_Halo_Boost_Interp = Use_Halo_Boost_Interp,
            Use_Halo_Boost_Extrap = Use_Halo_Boost_Extrap,
            Use_mdot_Interp = Use_mdot_Interp,
            Use_LowM_Extrap_mdot = Use_LowM_Extrap_mdot,
            Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot,
            Break_OmB_Lost_Ratio = OmB_Lost_Ratio * 1.1,
            Fix_mdot = Fix_mdot,
            Fixed_mdot = Fixed_mdot)
        lost_ratio = np.max(1 - r['omb_ratio'])
        dif = lost_ratio - OmB_Lost_Ratio

        return dif
    
    # For small bh (accretion inefficient) Solve might fail
    try:
        lg_fbh0_max = Solve(
            F = omb_dif,
            Xmin = lf_min,
            Xmax = 0,
            Precision = Precision,
            show_status = show_status)
    except:
        lg_fbh0_max = 0
    
    r = 10**lg_fbh0_max

    return r

#--------Part II: Curvature Perturbation----
# The part that actually deals with PBH formation, copied from my NanoGrav codes
try:
    from src.PBH_Formation import *
except:
    # for linked import
    from PBH_Formation import *

def Find_PsMax(
        kc = m2k(100),
        Sigma = 0.4,
        OmB_Lost_Ratio = 0.3,
        Use_EoR = True,
        DeltaC = 0.45,
        map_nx = 500,
        OmB0 = cosmo['OmB'],
        zmin = 8,
        zmax = 1000,
        nz = 5000,
        nm = 500,
        Use_Halo_Boost = False,
        sbh_range = [-3, 5],
        show_status = False,
        Use_Halo_Boost_Interp = True,
        Use_Halo_Boost_Extrap = True,
        Use_mdot_Interp = True,
        Use_LowM_Extrap_mdot = True,
        Use_HighM_Extrap_mdot = True,
        Precision = 1e-2,
        Show_Extrap_MSG = False,
        Fix_mdot = False,
        Fixed_mdot = 10
        ):
    '''
    Find maximum Ps amplitude
    OmB_Lost_Ratio : 
        maximumly allowed OmB lost ratio, default value of 0.3 is motivated by DES Low-Z results
        DOI: 10.1103/PhysRevD.105.023520
    Speed : 318s for nz = 1500, nm = 500, map_nx = 500
    '''
    
    # only consider Ps which give fbh in this range
    fbh_range = [1e-30, 0.95]

    def fbh_dif(LgA):
        # Find whether PBH produced by LgA over-accrets baryon
        '''
        MF = Phi_Profile(
            A = 10**LgA,
            Sigma = Sigma,
            kc = kc,
            DeltaC = DeltaC,
            map_nx = map_nx
        )
        
        # fbh above can be very big or small

        # Now find omb
        Evo = Get_Evo(
            mc = MF['mc'],
            sbh = MF['Width'],
            fbh0 = MF['fbh'],
            OmB0 = OmB0,
            Use_EoR = Use_EoR,
            zmin = zmin,
            zmax = zmax,
            nz = nz,
            Use_Halo_Boost = Use_Halo_Boost,
            nm = nm,
            sbh_range = sbh_range,
            show_status = False,
            LowM_method = LowM_method
        )
        '''
        Evo = Get_Evo(
            Use_Ps_domain = True,
            PsA = 10**LgA,
            kc = kc,
            SigmaPs = Sigma,
            DeltaC = DeltaC,
            OmB0 = OmB0,
            Use_EoR = Use_EoR,
            Use_Halo_Boost = Use_Halo_Boost,
            zmin = zmin,
            zmax = zmax,
            nz = nz,
            nm = nm,
            sbh_range = sbh_range,
            show_status = False,
            Use_Halo_Boost_Interp = Use_Halo_Boost_Interp,
            Use_Halo_Boost_Extrap = Use_Halo_Boost_Extrap,
            Use_mdot_Interp = Use_mdot_Interp,
            Use_LowM_Extrap_mdot = Use_LowM_Extrap_mdot,
            Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot,
            Break_OmB_Lost_Ratio = OmB_Lost_Ratio,
            Show_Extrap_MSG = Show_Extrap_MSG,
            Fix_mdot = Fix_mdot,
            Fixed_mdot = Fixed_mdot)
        
        omb_ratio = 1 - np.min(Evo['omb_ratio'])
        dif = omb_ratio - OmB_Lost_Ratio
        return dif

    # Smarter way to find scan region for A
    
    def Find_fbh_for_A(A):
        MF = Phi_Profile(
            A = A,
            Sigma = Sigma,
            kc = kc,
            DeltaC = DeltaC,
            map_nx = map_nx
        )
        return MF['fbh']

    def Find_A(fbh):
        # each call takes about 10s
        dif = lambda la: Find_fbh_for_A(10**la) - fbh
        try:
            r = Solve(F = dif, Xmin = -4, Xmax = 2, Precision = 1e-3)
        except:
            # Expand search
            r = Solve(F = dif, Xmin = -8, Xmax = 6, Precision = 1e-3)
        r = 10**r
        
        return r
    
    # Take 20s to optimise search range
    A_Search_Min = Find_A(fbh_range[0])
    A_Search_Max = Find_A(fbh_range[1])
    
    LA_Search_Min = np.log10(A_Search_Min)
    LA_Search_Max = np.log10(A_Search_Max)
    
    # It's possible for small PBHs to be unconstrained, if so return negative value
    try:
        LgA_max = Solve(F = fbh_dif, Xmin = LA_Search_Min, Xmax = LA_Search_Max, Precision = Precision, show_status = show_status)
        Amax = 10**LgA_max
        Pmax = Amax/(np.sqrt(2 * np.pi)*Sigma)
    except:
        warnings.warn('Solution not found, likely because BHs are too light, returning NaN, debug info:')
        print('A_Search_Min = ', A_Search_Min)
        print('A_Search_Max = ', A_Search_Max)
        print('OmB_dif @ Amin: ', fbh_dif(LA_Search_Min))
        print('OmB_dif @ Amax: ', fbh_dif(LA_Search_Max))

        Pmax = np.nan
    
    return Pmax
