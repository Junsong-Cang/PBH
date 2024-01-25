'''
Main module for computing PBH accretion history,
written by ZHZ & JSC
'''

from scipy import linalg, optimize
import numpy as np
from PyLab import *

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
    
def lambda_zzh(beta0,gamma0,exist_halo,r_B_eff,r_h,z,p, dlnx0 = -1e-4):
    #dlnx0 = -1e-4
    dlnx = dlnx0
    lnx_cut = 2.
    dlnx_cut = -1e-8

    x_h = r_h / r_B_eff
    if 1 > x_h:
        m_h_rBeff = 1.
    else:
        m_h_rBeff = (1 / x_h) ** p

    lglambda_max = 1.
    lglambda_min = -8.
    lglambda_h = (lglambda_max+lglambda_min)/2.
    #lambda_max = 10**lglambda_max
    #lambda_min = 10**lglambda_min
    lambda_h = 10**lglambda_h

    while 1:
        lnx_cut_change = 2.
        #boundary conditions
        lnx = 3.
        x = np.exp(lnx)

        rho = 1.
        #lnrho = np.log(rho)

        T = 1.
        lnT = np.log(T)

        v = lambda_h/rho/x**2
        lnv = np.log(v)

        while 1:
            if exist_halo == 1:#with DM halo
                if x > x_h:
                    m_h_r = 1.
                else:
                    m_h_r = (x/x_h)**p

                m_PBH = (1.+z)/3000.
                M_r = (m_PBH+m_h_r)/(m_PBH+m_h_rBeff)
            else:#no DM halo
                M_r = 1.

            A = v**2/T
            A0= (M_r/x-beta0 * v*x)/T

            B0= gamma0*x/v*(1./T-1)

            I = (A0 - B0 - 2 * A) / (5./3. - A)
            dlnv = - (2 - I) * dlnx
            #dlnrho = - I * dlnx
            dlnT = - (B0 + 2./3. * I) * dlnx
            #print('dlnv='+str(dlnv)+', dlnrho='+str(dlnrho)+', dlnT='+str(dlnT))
            lnx = lnx + dlnx
            lnv = lnv + dlnv
            #lnrho = lnrho + dlnrho
            lnT = lnT + dlnT
            x = np.exp(lnx)
            v = np.exp(lnv)
            #rho = np.exp(lnrho)
            T = np.exp(lnT)

            if lnx < lnx_cut_change:
                if dlnx>dlnx0 and dlnx<dlnx_cut:
                    dlnx = dlnx * 10**1
                lnx_cut_change = lnx_cut_change - 1

            if lnT < 2*lnv-0.510825623765991:# np.log(5./3.) = 0.510825623765991
                if lnx > lnx_cut and dlnx<dlnx_cut:
                    dlnx = dlnx*1e-1
                    break
                lglambda_max = lglambda_h
                #lambda_max = 10**lglambda_max
                dlnx = dlnx0
                break
            elif lnx < -8. or ( 10./3.*T/x - M_r/(x)**2 + gamma0 * (1-T)/v +beta0 * v )<0:
                if lnx>lnx_cut:# and dlnx < dlnx_cut:
                    dlnx = dlnx*1e-1
                    break
                lglambda_min = lglambda_h
                #lambda_min = 10**lglambda_min
                dlnx = dlnx0
                break

        lglambda_h = (lglambda_max+lglambda_min)/2
        lambda_h = 10**lglambda_h
        if lglambda_max - lglambda_min <0.01:
            break

    return 10**lglambda_max


def mdot_kernel(
        z = 10,
        M_PBH = 100,
        x_e = 1e-4,
        T_k = 10,
        Omega_b = 0.04897468161,
        Use_halo = 1,
        Use_feedback = 0,
        Use_EoR = 1):
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
    h = 0.6766
    Omb_LCDM = 0.04897468161
    k_B=1.3806542e-16
    m_H=1.66e-24 # H mass in g
    G=6.67259e-8 # cm^3/g/s^2
    M_SUN=1.989e33 # solar mass in g
    pc=3.086e18 # parsec in cm
    c=2.99792458e10
    sigma_T=6.65e-25
    m_e=9.1093897e-28 # electron mass in g
    rhoc=(2.7752e11)*h**2#M_SUN/Mpc^3
    Y_He=0.245
    n_H_z0=rhoc*Omb_LCDM*M_SUN*(1-Y_He)/m_H/(1e6*pc)**3
    n_He_z0=rhoc*Omb_LCDM*M_SUN*Y_He/(4*m_H)/(1e6*pc)**3

    p = 0.75 # power-law DM halo density profile
    rho_crit = 1.879e-29 * h**2 # g/cm^3
    rho_b = Omega_b * rho_crit * (1+z)**3 # g/cm^3
    L_Edd = 1.26e38*M_PBH # erg/s

    T_k_z = T_k
    x_e_z = x_e
    if Use_feedback == 0:
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

    #chi0 = r_Bh/r_h

    # lambda
    if Use_halo == 0:
        # naked PBH
        rho_CMB = 4.15e-5 * h**-2 * rho_crit * (1+z)**4 * c**2
        t_B = r_B * pc / (v_eff * 1e5)
        beta0  = 4 * x_e_z * sigma_T * rho_CMB * t_B / (3*m_H*c)
        gamma0 = 8 * x_e_z * sigma_T * rho_CMB * t_B / (3*m_e*c*(1+x_e_z))

        lambda_ad = 1./4.*(3./5.)**1.5
        lambda_iso = 1./4.*np.exp(3./2.)
        lambda0 = (lambda_ad+(lambda_iso-lambda_ad)*(gamma0**2/(88.+gamma0**2))**0.22 * np.exp(9./2./(3.+beta0**0.75)) / (np.sqrt(1.+beta0)+1.)**2 )/lambda_iso
        # lambda1 = lambda_zzh(beta0, gamma0, 0, r_B, 0, z, p)
        M_acc = 4 * np.pi * lambda0 * rho_b * (r_B*pc)**2 * (v_eff*1e5)#g/s

    elif Use_halo == 1:
        # PBH with DM halo
        r1 = r_B_eff_fun(z, M_PBH, v_eff, r_B, p)
        r_B_eff = r1[0]
        r_Bh = r1[1]
        r_h = r1[2]
        rho_CMB = 4.15e-5 * h ** -2 * rho_crit * (1 + z) ** 4 * c ** 2
        t_B_eff = r_B_eff * pc / (v_eff * 1e5)
        beta0_h = 4 * x_e_z * sigma_T * rho_CMB * t_B_eff / (3*m_H*c)
        gamma0_h = 8 * x_e_z * sigma_T * rho_CMB * t_B_eff / (3 * m_e * c * (1 + x_e_z))
        lambda1 = lambda_zzh(beta0_h, gamma0_h, 1, r_B_eff, r_h, z, p)
        M_acc = 4 * np.pi * lambda1 * rho_b * (r_B_eff*pc)**2 * (v_eff*1e5)#g/s

    # m_acc
    m_acc = M_acc * c**2 / L_Edd
    #m_acc = min(m_acc,10)
    return m_acc#, lambda1
