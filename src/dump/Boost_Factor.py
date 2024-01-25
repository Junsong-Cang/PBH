from p21c_tools import *
from src_1 import *

def Boost_Factor_Lite(z = 10,
                 OmB = cosmo['OmB'],
                 Use_EoR = False,
                 nmh = 1000,
                 show_status = 0,
                 ncpu = 1
                 ):
    '''
    Get fcoll and hmf integration
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
    r = Total_Boost

    return r

def BoostFactor(
        z = 10,
        m = 0.1,
        OmB = cosmo['OmB'],
        Use_EoR = False,
        nmh = 1000,
        Use_Halo_Profile_Interp = 0,
        show_status = 0,
        ncpu = 1,
        Use_Edd_Limit = True
    ):

    Is_Scalar(x = z, ErrorMethod = 2)
    OmC = cosmo['OmC']
    h = cosmo['h']
    OmM = OmB + OmC
    RhoCr = 2.775e11 * h**2 # msun/Mpc^3
    RhoM_cmv = RhoCr*OmM
    RhoC_cmv = RhoCr*OmC
    RhoB_avg = RhoCr*OmB*(1+z)**3
    
    mdot_ideal = Get_mdot(z = z, m = m, OmB = OmB, Use_EoR = Use_EoR, Use_Edd_Limit = False)
    
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
            Use_Interp = Use_Halo_Profile_Interp,
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
    
    return result
