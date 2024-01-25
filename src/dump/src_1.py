'''
Main module for computing PBH accretion history,
written by ZHZ & JSC
'''

from scipy import linalg
from sympy import *
from scipy import optimize
import numpy as np
import h5py, os, time, warnings, sys, h5py, platform
from PyLab import *
from Useful_Numbers import Cosmology as cosmo
from Useful_Numbers import Constants
from p21c_tools import *

if platform.system() == 'Darwin':
    main_path = '/Users/cangtao/cloud/GitHub/PBH/'

def r_B_eff_fun(z, M_PBH, v_eff, r_B, p):
    '''
    ---- inputs ----
    z : redshift
    M_PBH : PBH mass in msun
    v_eff: the velocity of the relative motion between the black hole and the ambient gas
    r_B: PBH Bondi radius
    p: power-law DM halo density profile
    '''
    G=6.67259e-8 # cm^3/g/s^2
    M_SUN=1.989e33 # solar mass in g
    pc=3.086e18 # parsec in cm

    M_h = 3000. * M_PBH / (1.+z) # DM halo mass
    r_h = 0.339 * 58. / (1.+z) * M_h**(1./3.) # DM halo radius
    r_Bh = G * M_SUN * M_h / (v_eff * 1e5)**2 / pc # DM halo's Bondi radius
    r_B_eff = r_Bh + r_B

    if r_B_eff < r_h:
        def f(x):
            return r_B/x + r_Bh/x*(x/r_h)**p - 1
        r_B_eff=optimize.brentq(f, r_B, r_h, rtol=1e-4)
    r1 = np.array([r_B_eff, r_Bh, r_h])
    return r1

def lambda_zzh(beta0,gamma0,exist_halo,r_B_eff,r_h,z,p):
    dlnx = -1e-4
    lglambda_max = 4.
    lglambda_min = -8.
    lglambda_h = (lglambda_max+lglambda_min)/2.
    lambda_max = 10**lglambda_max
    lambda_min = 10**lglambda_min
    lambda_h = 10**lglambda_h

    while 1:
        #boundary conditions
        lnx = 3.
        x = np.exp(lnx)

        rho = 1.
        lnrho = np.log(rho)

        T = 1.
        lnT = np.log(T)

        v = lambda_h/rho/x**2
        lnv = np.log(v)

        while 1:
            if exist_halo == 1:  # with DM halo
                x_h = r_h / r_B_eff
                lnx_h = np.log(x_h)
                if lnx > lnx_h:
                    m_h_r = 1.
                else:
                    m_h_r = (np.exp(lnx - lnx_h)) ** p
                if 0 > lnx_h:
                    m_h_rBeff = 1.
                else:
                    m_h_rBeff = (np.exp(-lnx_h)) ** p
                m_PBH = (1. + z) / 3000.
                M_r = (m_PBH + m_h_r) / (m_PBH + m_h_rBeff)
            else:  # no DM halo
                M_r = 1.

            if gamma0 < 10:
                A1 = np.exp(2 * lnv)
                A2 = np.exp(lnT)
                A3 = np.exp(lnT)
                A4 = M_r / np.exp(lnx) - beta0 * np.exp(lnv + lnx)
                B1 = 0.
                B2 = 2. / 3. * np.exp(lnv)
                B3 = -np.exp(lnv)
                B4 = -gamma0 * np.exp(lnx) * (np.exp(-lnT) - 1)
                C1 = 1.
                C2 = 1.
                C3 = 0.
                C4 = 2.

                I1 = np.array([[A1, A2, A3],
                               [B1, B2, B3],
                               [C1, C2, C3]])
                I2 = np.array([A4, B4, C4]) * (-dlnx)
                X = linalg.solve(I1, I2)

                dlnv = X[0]
                dlnrho = X[1]
                dlnT = X[2]

                lnx = lnx + dlnx
                lnv = lnv + dlnv
                lnrho = lnrho + dlnrho
                lnT = lnT + dlnT
            else:
                A1 = np.exp(2 * lnv)
                A2 = np.exp(lnT)
                A4 = M_r / np.exp(lnx) - beta0 * np.exp(lnv + lnx)
                C1 = 1.
                C2 = 1.
                C4 = 2.

                I1 = np.array([[A1, A2],
                               [C1, C2]])
                I2 = np.array([A4, C4]) * (-dlnx)
                X = linalg.solve(I1, I2)

                dlnv = X[0]
                dlnrho = X[1]
                dlnT = 0

                lnx = lnx + dlnx
                lnv = lnv + dlnv
                lnrho = lnrho + dlnrho
                lnT = lnT + dlnT


            if -np.exp(lnv)-5./3.*np.exp(lnT)/(-np.exp(lnv)) < 0:
                lglambda_max = lglambda_h
                lambda_max = 10**lglambda_max
                break
            elif lnx < -8. or ( 10./3.*np.exp(lnT)/np.exp(lnx) - M_r/(np.exp(lnx))**2 + gamma0 * (1-np.exp(lnT))/np.exp(lnv) +beta0 * np.exp(lnv) )<0:
                lglambda_min = lglambda_h
                lambda_min = 10**lglambda_min
                break

        lglambda_h = (lglambda_max+lglambda_min)/2
        lambda_h = 10**lglambda_h
        if lglambda_max - lglambda_min <0.01:
            break

    return lambda_max

def mdot(z = 10,
         M_PBH = 100,
         x_e = 1e-4,
         T_k = 10,
         Omega_b = 0.0484,
         Use_halo = 1,
         Use_EoR = 1,
         Use_feedback = 0):
    '''
    Dimensionless accretion rate for PBH
    ---- inputs ----
    z : redshift
    M_PBH : PBH mass in msun
    x_e : ionisation fraction at z
    T_k : gas temterature at z
    Omega_b : baryon density fraction
    Use_halo : whether or not to use dm halo
    Use_feedback : whether or not to use feedback in x_e and T_k,
                   if False then x_e and T_k inputs are ignored and LCDM is used
    '''
    # parameters
    h = cosmo['h']
    k_B=1.3806542e-16
    m_H=1.66e-24 # H mass in g
    G=6.67259e-8 # cm^3/g/s^2
    M_SUN=1.989e33 # solar mass in g
    pc=3.086e18 # parsec in cm
    c=2.99792458e10
    sigma_T=6.65e-25
    m_e=9.1093897e-28 # electron mass in g
    rhoc=(2.7752e11)*h**2#M_SUN/Mpc^3
    Y_He=0.24
    n_H_z0=rhoc*Omega_b*M_SUN*(1-Y_He)/m_H/(1e6*pc)**3
    n_He_z0=rhoc*Omega_b*M_SUN*Y_He/(4*m_H)/(1e6*pc)**3
    
    p = 0.75 # power-law DM halo density profile
    rho_crit = 1.879e-29 * h**2 # g/cm^3
    rho_b = Omega_b * rho_crit * (1+z)**3 # g/cm^3
    L_Edd = 1.26e38*M_PBH # erg/s

    T_k_z = T_k
    x_e_z = x_e
    if not Use_feedback:
        x_e_z, T_k_z = LCDM_HyRec(z = z, Use_EoR = Use_EoR)

    n_H = n_H_z0 * (1 + z) ** 3
    n_He = n_He_z0 * (1 + z) ** 3
    n_tot = n_H * (1 + x_e_z) + n_He
    mu = rhoc * Omega_b * M_SUN / (1e6 * pc) ** 3 * (1 + z) ** 3 / n_tot / m_H

    cs = np.sqrt( 5./3.*k_B*T_k_z/(mu*m_H))/1e5 # km/s
    v_L = np.min([1., z / 1e3]) * 30.#km/s
    v_eff = np.max([cs,np.sqrt(cs*v_L)])#km/s

    # Bondi radius
    r_B = G * M_SUN * M_PBH / (v_eff*1e5)**2/pc#pc
    r1 = r_B_eff_fun(z,M_PBH,v_eff,r_B,p)
    r_B_eff = r1[0]
    r_Bh = r1[1]
    r_h = r1[2]

    # lambda
    if Use_halo == 0:
        # naked PBH
        rho_CMB = 4.15e-5 * h**-2 * rho_crit * (1+z)**4 * c**2
        t_B = r_B * pc / (v_eff * 1e5)
        beta0  = 4 * x_e_z * sigma_T * rho_CMB * t_B / (3*m_H*c)
        gamma0 = 8 * x_e_z * sigma_T * rho_CMB * t_B / (3*m_e*c*(1+x_e_z))

        lambda_ad = 1./4.*(3./5.)**1.5
        lambda_iso = 1./4.*np.exp(3./2.)
        lambda0 = (lambda_ad+(lambda_iso-lambda_ad)*(gamma0**2/(88.+gamma0**2))**0.22 * exp(9./2./(3.+beta0**0.75)) / (np.sqrt(1.+beta0)+1.)**2 )/lambda_iso
        M_acc = 4 * np.pi * lambda0 * rho_b * (r_B*pc)**2 * (v_eff*1e5)#g/s

    elif Use_halo == 1:
        # PBH with DM halo
        rho_CMB = 4.15e-5 * h ** -2 * rho_crit * (1 + z) ** 4 * c ** 2
        t_B_eff = r_B_eff * pc / (v_eff * 1e5)
        beta0_h = 4 * x_e_z * sigma_T * rho_CMB * t_B_eff / (3*m_H*c)
        gamma0_h = 8 * x_e_z * sigma_T * rho_CMB * t_B_eff / (3 * m_e * c * (1 + x_e_z))
        lambda_h = lambda_zzh(beta0_h, gamma0_h, 1, r_B_eff, r_h, z, p)
        M_acc = 4 * np.pi * lambda_h * rho_b * (r_B_eff*pc)**2 * (v_eff*1e5)#g/s

    # m_acc
    m_acc = M_acc * c**2 / L_Edd
    # do this at Get_mdot
    #if(m_acc >10.0):
    #   m_acc=10.0

    return m_acc

def Get_BH_Interp_Data():
    '''
    Get interpolation data
    data produced by Get_Tables.py and Get_Tables_Booster.py
    '''
    if platform.system() == 'Darwin':
        # running on mac
        DataPath = '/Users/cangtao/cloud/GitHub/PBH/data/'
    else:
        # running on IHEP server
        DataPath = '/home/dm/watson/work/Accreting_PBH/data/'
    
    # First check HighZ file, useful for recombination calculation
    # Contains fields: mdot_vec, z, m
    # mdot_vec index: [mass, z]

    HighZ_File = DataPath + 'mdot_table_HighZ.npz'
    HighZ_Tab = np.load(HighZ_File)

    def Load_LCDM_Table():
    
        # Small BH
        TabFile = DataPath + 'mdot_table.h5'
        f = h5py.File(TabFile, 'r')
        mdot_table = f['mdot_table'][:]
        m_axis = f['m_axis'][:]
        z_axis = f['z_axis'][:]
        f.close()
        Light_BH_Table = {'mdot_table' : mdot_table,
                          'm_axis' : m_axis,
                          'z_axis' : z_axis}
    
        # Even smaller BH
        TabFile = DataPath + 'mdot_table_tiny.h5'
        f = h5py.File(TabFile, 'r')
        mdot_tiny = f['mdot_table'][:]
        m_tiny = f['m_axis'][:]
        z_tiny = f['z_axis'][:]
        f.close()
        Tiny_BH_Table = {
            'mdot_table' : mdot_tiny,
            'm_axis' : m_tiny,
            'z_axis' : z_tiny}

        # Heavy BH
        TabFile = DataPath + 'mdot_table_booster.h5'
        f = h5py.File(TabFile, 'r')
        mdot_table_booster = f['mdot_table'][:] # index: [zid, mid]
        Booster_Shape = np.shape(mdot_table_booster)
        nm = Booster_Shape[1]
        z_vec = f['z_vec'][:]
        mmin = np.array( f['mmin'] )
        mmax = f['mmax'][:]
        f.close()
    
        Booster_BH_Table = {
            'mdot_table': mdot_table_booster,
            'z_vec' : z_vec,
            'mmin' : mmin,
            'mmax' : mmax,
            'nm' : nm}
        r = {
            'Tiny': Tiny_BH_Table,
            'Light': Light_BH_Table,
            'Booster': Booster_BH_Table
            }

        return r
    
    def Load_EoR_Table():
    
        # Small BH
        TabFile = DataPath + 'mdot_table_EoR.h5'
        f = h5py.File(TabFile, 'r')
        mdot_table = f['mdot_table'][:]
        m_axis = f['m_axis'][:]
        z_axis = f['z_axis'][:]
        f.close()
        Light_BH_Table = {'mdot_table' : mdot_table,
                          'm_axis' : m_axis,
                          'z_axis' : z_axis}
    
        # Even smaller BH
        TabFile = DataPath + 'mdot_table_tiny_EoR.h5'
        f = h5py.File(TabFile, 'r')
        mdot_tiny = f['mdot_table'][:]
        m_tiny = f['m_axis'][:]
        z_tiny = f['z_axis'][:]
        f.close()
        Tiny_BH_Table = {
            'mdot_table' : mdot_tiny,
            'm_axis' : m_tiny,
            'z_axis' : z_tiny}

        # Heavy BH
        TabFile = DataPath + 'mdot_table_booster_EoR.h5'
        f = h5py.File(TabFile, 'r')
        mdot_table_booster = f['mdot_table'][:] # index: [zid, mid]
        Booster_Shape = np.shape(mdot_table_booster)
        nm = Booster_Shape[1]
        z_vec = f['z_vec'][:]
        mmin = np.array( f['mmin'] )
        mmax = f['mmax'][:]
        f.close()
    
        Booster_BH_Table = {
            'mdot_table': mdot_table_booster,
            'z_vec' : z_vec,
            'mmin' : mmin,
            'mmax' : mmax,
            'nm' : nm}
        
        r = {
            'Tiny': Tiny_BH_Table,
            'Light': Light_BH_Table,
            'Booster': Booster_BH_Table
            }

        return r
    
    # Boost Factor Table, author: Boost_Tab.py
    datafile = DataPath + 'Boost_Tab.npz'
    Boost_Factor_Tab = np.load(datafile)

    # Put them all in here
    LCDM_Tab = Load_LCDM_Table()
    EoR_Tab = Load_EoR_Table()
    r = {
        'HighZ_Tab' : HighZ_Tab,
        'LCDM_Tab' : LCDM_Tab, 
        'EoR_Tab' : EoR_Tab,
        'Boost_Factor_Tab': Boost_Factor_Tab}
    return r

# Use as global variable for speed (avoid reload)
BH_mdot_interp_table = Get_BH_Interp_Data() 

def Get_mdot(
        z = 10,
        m = 100,
        OmB = cosmo['OmB'],
        Do_not_Interpolate = False,
        Use_LowM_Extrap = True,
        Use_HighM_Extrap = True,
        Show_Extrap_Warning = True,
        Use_EoR = True,
        Use_Edd_Limit = True):
    '''
    A fast interface for mdot,
    use interpolation for low masses to speed up
    Interpolation precision:
    0.5% on average, maximum error ~ 8% (for [nm, nz] = [100, 50])
    Interpolation speed:
    20000 calls per second (1 cpu, MBP M2 Ultra)
    ---- inputs ----
    z : redshift
    m : bh mass in msun
    Do_not_Interpolate : abort interpolation, useful for testing
    Use_LowM_Extrap : use extrapolation for LowM overflow
    Use_HighM_Extrap : use extrapolation for HighM overflow
    Show_Extrap_Warning : Complaing when extrapolation is used
    '''
    
    Small = 1e-250
    EoR_21cmFAST_MaxZ = 36
    OmB_LCDM = cosmo['OmB']
    
    # Allowed overflow fraction, will be useful for later validation
    Extrap_Redundancy = 1e5
    
    # Set a floor for OmB_Ratio, sometimes in numerical integration OmB may transite from Small to -Small and this will lead to errors
    OmB_Ratio = OmB/OmB_LCDM
    if OmB_Ratio < Small:
        return 0

    if Use_EoR and z < EoR_21cmFAST_MaxZ:
        mdot_table = BH_mdot_interp_table['EoR_Tab']
    else:
        mdot_table = BH_mdot_interp_table['LCDM_Tab']
    
    Tiny_BH_Table = mdot_table['Tiny']
    Light_BH_Table = mdot_table['Light']
    Booster_BH_Table = mdot_table['Booster']
    HighZ_Tab = BH_mdot_interp_table['HighZ_Tab']

    m_axis_lite = Light_BH_Table['m_axis']
    mdot_tab_lite = Light_BH_Table['mdot_table']
    zp_axis_lite = Light_BH_Table['z_axis'] + 1

    m_axis_tiny = Tiny_BH_Table['m_axis']
    mdot_tab_tiny = Tiny_BH_Table['mdot_table']
    zp_axis_tiny = Tiny_BH_Table['z_axis'] + 1

    z_axis_booster = Booster_BH_Table['z_vec']

    m_axis_HZ = HighZ_Tab['m']
    mdot_tab_HZ = HighZ_Tab['mdot_vec']
    zp_axis_HZ = HighZ_Tab['z'] + 1

    zp = z + 1
    
    # Initialise mbh
    mbh = m
    Redundancy = 1e-5

    # Check interpolation table
    
    # ----HighZ----
    # For z = [1000, 2000], m = [1e-3, 1e5]
    if (not Do_not_Interpolate) and Within_Range(zp, zp_axis_HZ):
        
        if not Within_Range(mbh, m_axis_HZ):
            if mbh < m_axis_HZ[0]:
                if Use_LowM_Extrap:
                    mbh_ = m_axis_HZ[0] * (1+Redundancy)
                    OverFlow = mbh_/m
                    if OverFlow < Extrap_Redundancy:
                        mbh = mbh_
                        if Show_Extrap_Warning:
                            #warnings.warn('Use_LowM_Extrap triggered, info:')
                            print('Use_LowM_Extrap triggered at HighZ module, input : [z, m] = [', '{0:.3f}'.format(z), ',',  '{0:.4e}'.format(m), ']',  
                                    ', interp MinM: ', '{0:.4e}'.format(m_axis_HZ[0]), ', Overflow :', '{0:.4e}'.format(OverFlow))
            else:
                if Use_HighM_Extrap:
                    mbh_ = m_axis_HZ[-1] * (1-Redundancy)
                    OverFlow = m/mbh_
                    if OverFlow < Extrap_Redundancy:
                        mbh = mbh_
                        if Show_Extrap_Warning:
                            # warnings.warn('Use_HighM_Extrap triggered, info:')
                            print('Use_HighM_Extrap triggered at HighZ module, input : [z, m] = [', '{0:.3f}'.format(z), ',',  '{0:.4e}'.format(m), ']',  
                                    ', interp MaxM: ', '{0:.4e}'.format(m_axis_HZ[-1]), ', Overflow :', '{0:.4e}'.format(OverFlow))

        if Within_Range(mbh, m_axis_HZ):
            r = Interp_2D(
                Tab = mdot_tab_HZ,
                x_axis = m_axis_HZ,
                y_axis = zp_axis_HZ,
                x_target = mbh,
                y_target = zp,
                Use_Log_X = True,
                Use_Log_Y = True,
                Use_Log_Z = True
            )

            if Use_Edd_Limit:
                r = min(OmB_Ratio * r, 10)
            else:
                r = OmB_Ratio * r
            return r

    # ----Tiny----
    # For m = [1e-3, 1e0]
    if Within_Range(zp, zp_axis_tiny) and not Do_not_Interpolate:
        
        # HighM overflow will be diverted to Lite or Booster Tables
        if mbh < m_axis_tiny[0] and Use_LowM_Extrap:
            mbh_ = m_axis_tiny[0] * (1+Redundancy)
            OverFlow = mbh_/m
            if OverFlow < Extrap_Redundancy:
                mbh = mbh_
                if Show_Extrap_Warning:
                    warnings.warn('Use_LowM_Extrap triggered, info:')
                    print('Use_LowM_Extrap triggered at Tiny module, input : [z, m] = [', '{0:.3f}'.format(z), ',',  '{0:.4e}'.format(m), ']',  
                        ', interp MinM: ', '{0:.4e}'.format(m_axis_tiny[0]), ', Overflow :', '{0:.4e}'.format(OverFlow))

        if Within_Range(mbh, m_axis_tiny):
            r = Interp_2D(
                Tab = mdot_tab_tiny, 
                x_axis = m_axis_tiny, 
                y_axis = zp_axis_tiny,
                x_target = mbh,
                y_target = zp,
                Use_Log_X = True,
                Use_Log_Y = True,
                Use_Log_Z = True,
                )
            if Use_Edd_Limit:
                r = min(OmB_Ratio * r, 10)
            else:
                r = OmB_Ratio * r
            return r

    # ----LowM----
    # For m = [1e0, 1e4]
    # Extrap not required as overflow will be diverted to Booster Table

    if (Within_Range(mbh, m_axis_lite) and Within_Range(zp, zp_axis_lite)) and not Do_not_Interpolate:

        r = Interp_2D(
            Tab = mdot_tab_lite, 
            x_axis = m_axis_lite, 
            y_axis = zp_axis_lite,
            x_target = mbh,
            y_target = zp,
            Use_Log_X = True,
            Use_Log_Y = True,
            Use_Log_Z = True,
            )
        if Use_Edd_Limit:
            r = min(OmB_Ratio * r, 10)
        else:
            r = OmB_Ratio * r
        return r

    # ----Booster----
    # For m > 1e4
    
    if Within_Range(z, z_axis_booster) and not Do_not_Interpolate:
        zid1 = Find_Index(z, z_axis_booster)
        zid2 = zid1 + 1
        
        mdot_vec_1 = Booster_BH_Table['mdot_table'][zid1,:]
        mdot_vec_2 = Booster_BH_Table['mdot_table'][zid2,:]
        
        # Get m axis
        mmax_1 = Booster_BH_Table['mmax'][zid1]
        mmax_2 = Booster_BH_Table['mmax'][zid2]
        Lg_mmax_1 = np.log10(mmax_1)
        Lg_mmax_2 = np.log10(mmax_2)
        lg_mmin = np.log10(Booster_BH_Table['mmin'])

        m_vec_1 = np.logspace(lg_mmin, Lg_mmax_1, Booster_BH_Table['nm'])
        m_vec_2 = np.logspace(lg_mmin, Lg_mmax_2, Booster_BH_Table['nm'])
        lm_vec_1 = np.log10(m_vec_1)
        lm_vec_2 = np.log10(m_vec_2)

        # donno whether mbh is in range
        if (not Within_Range(mbh, m_vec_1)) or (not Within_Range(mbh, m_vec_2)):
            # This can happen at HighM or at very very low M
            if Use_HighM_Extrap and (mbh > m_vec_1[0]):
                MaxM_1 = m_vec_1[-1]
                MaxM_2 = m_vec_2[-1]
                MinM = min(MaxM_1, MaxM_2)
                mbh_ = MinM * (1-Redundancy)

                OverFlow = m/mbh_
                if OverFlow < Extrap_Redundancy:
                    mbh = mbh_
                    if Show_Extrap_Warning:
                        warnings.warn('Use_HighM_Extrap triggered, info:')
                        print('Use_HighM_Extrap triggered at Booster module, input : [z, m] = [', '{0:.3f}'.format(z), ',',  '{0:.4e}'.format(m), ']',  
                              ', interp MaxM: ', '{0:.4e}'.format(MinM), ', Overflow :', '{0:.4e}'.format(OverFlow))
        
        if Within_Range(mbh, m_vec_1) and Within_Range(mbh, m_vec_2):
            LgM = np.log10(mbh)
            # find mdot for z1 and z2
            mdot_z1 = np.interp(x = LgM, xp = lm_vec_1, fp = mdot_vec_1)
            mdot_z2 = np.interp(x = LgM, xp = lm_vec_2, fp = mdot_vec_2)
            # interpolate z axis
            x1 = np.log10(z_axis_booster[zid1] + 1)
            x2 = np.log10(z_axis_booster[zid2] + 1)
            x = np.log10(zp)

            r = (mdot_z2 - mdot_z1)*(x - x1)/(x2 - x1) + mdot_z1
            if Use_Edd_Limit:
                r = min(OmB_Ratio * r, 10)
            else:
                r = OmB_Ratio * r
            return r
    
    if not Do_not_Interpolate:
        if Use_HighM_Extrap and Use_LowM_Extrap:
            warnings.warn('Extrap failed to catch overflow within allowed range, I am handling this to mdot(). Overflow:')
            # print(OverFlow)
        else:
            warnings.warn('No match found in interpolation table, diverting to mdot(). Possible cause: Extrapolation deactivated, severe overflow')

    # If not in range and not Extrap then following lines will be executed
    # the inputs bellow must conform with interpolation table
    r = mdot(
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
        Use_Halo_Profile_Interp = 0,
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
        Boost_Table = BH_mdot_interp_table['Boost_Factor_Tab']
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
        Use_HighM_Extrap_mdot = True
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
    '''

    msun = Constants['msun']
    c = Constants['c']

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
            Use_HighM_Extrap = Use_HighM_Extrap_mdot
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
    
    H = Hubble(z)
    r = dmdt/H, Lx
        
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
        nz = 10000,
        Use_Halo_Boost = False,
        Use_EoR = True,
        Use_OmB_FeedBack = True,
        Fix_mdot = False,
        Fixed_mdot = 10,
        Use_Halo_Boost_Interp = True,
        Use_Halo_Boost_Extrap = True,
        Use_mdot_Interp = True,
        Use_LowM_Extrap_mdot = True,
        Use_HighM_Extrap_mdot = True
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

    m_vec[0] = m0
    OmB_vec[0] = OmB0
    dmdx_vec[0] = 0
    Lbh_vec[0] = 0
    
    # Start evolving
    for xid in np.arange(1, nz):
        xid_prev = xid-1

        x_prev = x_vec[xid_prev]
        OmB_prev = OmB_vec[xid_prev]
        m_prev = m_vec[xid_prev]
        if OmB_prev < Small:
            dm_dx, Lbh = 0.0, 0.0
        else:
            dm_dx, Lbh = dmdx(
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
                Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot)
        
        dmdx_vec[xid] = dm_dx
        Lbh_vec[xid] = Lbh

        if Use_OmB_FeedBack and not Fix_mdot:
            dOmB_dx = dOmBdx(
                dm_dx = dm_dx,
                fbh0 = fbh0,
                m0 = m0
            )
        else:
            dOmB_dx = 0

        dm = dm_dx * dx
        dOmB = dOmB_dx * dx
        
        m_vec[xid] = m_prev + dm
        OmB_vec[xid] = max(OmB_prev + dOmB, Small / 2)
    
    # Get dEdVdt
    OmC = cosmo['OmC']
    h = cosmo['h']
    RhoCr = 2.775e11 * h**2
    RhoC = RhoCr * OmC *(1+z_vec)**3 # mcun/Mpc^3
    nbh = fbh0 * RhoC/m0 # Mpc^3
    
    Q = Constants['Q']
    Mpc_cm = Constants['Mpc']*100
    
    Lbh_vec_eV = Lbh_vec/Q
    dEdVdt = nbh * Lbh_vec
    # Convert from J/s/Mpc^3 to eV/s/cm^3:
    dEdVdt = dEdVdt / Q / Mpc_cm**3

    r = {'z' : z_vec, 'm_ratio' : m_vec/m0, 'OmB_ratio' : OmB_vec/OmB0, 'fbh_ratio' : m_vec/m0, 
        'dmdx_vec' : dmdx_vec, 'x_vec' : x_vec, 'dEdVdt' : dEdVdt, 'Lbh_vec_eV' : Lbh_vec_eV}
    return r

def Find_fbhmax_lite(
        m0 = 1e3,
        OmB0 = cosmo['OmB'],
        OmB_Lost_Ratio = 0.2,
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
        Use_HighM_Extrap_mdot = True
        ):
    '''
    Find fbh limit for monochromatic BHs
    ----inputs----
    m0 : initial mass in msun
    OmB0 : Initial OmB
    OmB_Lost_Ratio : maximumly allowed OmB lost ratio
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
        dEdVdt_zmax = 1001,
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
        Fix_mdot = False,
        Fixed_mdot = 10
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

    '''
    
    OmC = cosmo['OmC']
    h = cosmo['h']
    Q = Constants['Q']
    Mpc_cm = Constants['Mpc']*100

    RhoCr = 2.775e11 * h**2 # msun/Mpc^3
    RhoC_cmv = RhoCr*OmC
    
    Small = 1e-280
    if dEdVdt_zmax < zmax:
        raise Exception('dEdVdt_zmax is smaller than zmax')
        
    # Get x array for time-evolution
    x1 = z2x(dEdVdt_zmax)
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
        
        if (f0 < 1e-60) or (np.isnan(SigmaBH)):
            warnings.warn('Input Ps generates no PBH (fbh<1e-60), I am returning NaN')
            return np.nan
        
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
    fbh_ratio = np.linspace(0, 1, nz)
    dEdVdt = np.linspace(0, 1, nz)
    
    # Initial Condition
    m_vec[0,:] = m0_vec[:]
    omb_vec[0] = OmB0
    fbh_ratio[0] = 1
    dEdVdt[0] = 0

    if show_status:
        print(' Mass range: [', m0_vec[0], ', ', m0_vec[-1])

    RhoBH0_cmv = f0 * RhoC_cmv
    Phi_vec[0,:] = Phi0[:]
    
    for zid in np.arange(1, nz):
        if show_status:
            print('status = ', zid/nz)
        zid_prev = zid-1
        x_prev = x_vec[zid_prev]
        omb_prev = omb_vec[zid_prev]
        OmB_Lost_Fraction = (OmB0 - omb_prev)/OmB0
        
        # Find mass evo
        for mid in np.arange(0, nm):
            try:
                if omb_prev < Small:
                    # Numerical cut, this means there is no baryon left
                    dm_dx, Lbh = 0.0, 0.0
                elif OmB_Lost_Fraction > 1.05*Break_OmB_Lost_Ratio:
                    # Economical cut
                    dm_dx, Lbh = 0.0, 0.0
                elif zid>1 and z_vec[zid] > zmax:
                    dm_dx = 0
                    Lbh =Lbh_vec[zid_prev,mid]
                else:
                    dm_dx, Lbh = dmdx(
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
                        Fixed_mdot = Fixed_mdot
                        )
            except:
                print('----crash iminent----')
                print('z = ', x2z(x_prev), ', m = ', m_vec[zid_prev, mid], ', OmB_ratio = ', omb_prev/OmB0)
                raise Exception('dmdx crashed')
            
            dm = dm_dx*dx
            m_vec[zid, mid] = m_vec[zid_prev, mid] + dm
            m_ratio[zid, mid] = m_vec[zid, mid]/m0_vec[mid]
            Lbh_vec[zid, mid] = Lbh
            

            # Find dEdVdt for mono----
            nbh = RhoBH0_cmv * ((1+z_vec[zid])**3)/m0_vec[mid] # Mpc^-3
            swap = nbh * Lbh # J/s/Mpc^3
            swap = swap/Q/Mpc_cm**3 # eV/s/cm^3
            dEdVdt_mono[mid] = swap

        dEdVdt[zid] = np.trapz(x = lnm0, y = dEdVdt_mono * Phi0)

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
        'Lbh_vec' : Lbh_vec
        }

    return r

def Find_fbhmax(
        mc = 1e2,
        sbh = 0.5,
        OmB_Lost_Ratio = 0.2,
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
        Use_HighM_Extrap_mdot = True
    ):
    '''
    Find maximum fbh0, takes about 14 Get_Evo calls
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
            Use_Halo_Boost_Interp = Use_Halo_Boost_Interp,
            Use_Halo_Boost_Extrap = Use_Halo_Boost_Extrap,
            Use_mdot_Interp = Use_mdot_Interp,
            Use_LowM_Extrap_mdot = Use_LowM_Extrap_mdot,
            Use_HighM_Extrap_mdot = Use_HighM_Extrap_mdot,
            Break_OmB_Lost_Ratio = OmB_Lost_Ratio * 1.1
        )
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
        OmB_Lost_Ratio = 0.1,
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
        Precision = 1e-2
        ):
    '''
    Find maximum Ps amplitude
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
            Use_Ps_domaim = True,
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
        )
        
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
        warnings.warn('Solution not found, likely because BHs are too light, returning NaN')
        Pmax = np.nan
    
    return Pmax
