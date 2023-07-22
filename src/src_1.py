'''
Main module for computing PBH accretion history,
written by ZHZ & JSC
'''

from scipy import linalg
from sympy import *
from scipy import optimize
from PyLab import *
from Useful_Numbers import *
import os, sys, h5py
from joblib import Parallel, delayed

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
    lglambda_max = 1.
    lglambda_min = -4.
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
            if exist_halo == 1:#with DM halo
                x_h = r_h/r_B_eff
                lnx_h = np.log(x_h)
                if lnx > lnx_h:
                    m_h_r = 1.
                else:
                    m_h_r = (np.exp(lnx-lnx_h))**p
                if 0 > lnx_h:
                    m_h_rBeff = 1.
                else:
                    m_h_rBeff = (np.exp(-lnx_h))**p
                m_PBH = (1.+z)/3000.
                M_r = (m_PBH+m_h_r)/(m_PBH+m_h_rBeff)
            else:#no DM halo
                M_r = 1.

            A1 = np.exp(2*lnv)
            A2 = np.exp(lnT)
            A3 = np.exp(lnT)
            A4 = M_r/np.exp(lnx)-beta0 * np.exp(lnv+lnx)
            B1 = 0.
            B2 = 2./3.*np.exp(lnv)
            B3 = -np.exp(lnv)
            B4 = -gamma0*np.exp(lnx)*(np.exp(-lnT)-1)
            C1 = 1.
            C2 = 1.
            C3 = 0.
            C4 = 2.

            I1 = np.array([[A1, A2, A3],
                           [B1, B2, B3],
                           [C1, C2, C3]])
            I2 = np.array([A4, B4, C4]) * (-dlnx)
            X = linalg.solve(I1,I2)

            dlnv = X[0]
            dlnrho = X[1]
            dlnT = X[2]

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
    h=0.6781
    Omb_LCDM=0.02226/h**2
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
    n_H_z0=rhoc*Omb_LCDM*M_SUN*(1-Y_He)/m_H/(1e6*pc)**3
    n_He_z0=rhoc*Omb_LCDM*M_SUN*Y_He/(4*m_H)/(1e6*pc)**3
    
    p = 0.75 # power-law DM halo density profile
    rho_crit = 1.879e-29 * h**2 # g/cm^3
    rho_b = Omega_b * rho_crit * (1+z)**3 # g/cm^3
    L_Edd = 1.26e38*M_PBH # erg/s

    T_k_z = T_k
    x_e_z = x_e
    if Use_feedback == 0:
        x_e_z, T_k_z = LCDM_HyRec(z)

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
        lambda0 = (lambda_ad+(lambda_iso-lambda_ad)*(gamma0**2/(88.+gamma0**2))**0.22 * exp(9./2./(3.+beta0**0.75)) / (np.sqrt(1.+beta0)+1.)**2 )/lambda_iso
        #lambda1 = lambda_zzh(beta0, gamma0, 0, r_B, 0, z, p)
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
    if(m_acc >10.0):
       m_acc=10.0

    return m_acc

'''
# Get interpolation table, written by JSC

reload = 0
ncpu = 12
TabFile = '/Users/cangtao/cloud/GitHub/PBH/data/mdot_table.h5'
m_vec = np.logspace(0, 4, 100)
zp_vec = np.logspace(np.log10(6), 3.2, 50)

def Get_mdot_vec(m):
    nz = len(zp_vec)
    r = np.linspace(0, 1, nz)
    # os.system('echo ---- >> /Users/cangtao/Desktop/tmp.txt')
    for idx in np.arange(0, nz):
        os.system('echo ---- >> /Users/cangtao/Desktop/tmp.txt')
        z = zp_vec[idx] - 1
        try:
            r[idx] = mdot(z = z, M_PBH = m, x_e = 1e-3, T_k = 100, Use_halo = 1, Use_feedback = 0)
        except:
            print('--------I have a bad feeling about this--------')
            print('r = ', r[idx])
            print('z = ', z)
            print('m = ', m)
            raise Exception('something went wrong ')
    return r

if reload:
    t1 = TimeNow()
    mdot_vec = Parallel(n_jobs=ncpu)(delayed(Get_mdot_vec)(x) for x in m_vec)
    # np.savez('tmp.npz', mdot_vec = mdot_vec)
    np.savez('tmp.npz', mdot_vec = mdot_vec)
    Timer(t1)

# saving data to h5 for faster loading
mdot_vec = np.load('tmp.npz')['mdot_vec']

f = h5py.File(TabFile, 'w')
f.create_dataset('mdot_table', data = mdot_vec)
f.create_dataset('m_axis', data = m_vec)
f.create_dataset('z_axis', data = zp_vec - 1)
f.close()

f = h5py.File(TabFile, 'r')
xx = f['mdot_table'][:]
f.close()

print(np.shape(xx))
'''

def Get_BH_Interp_Data():
    '''
    Get interpolation data for small BH
    '''
    TabFile = '/Users/cangtao/cloud/GitHub/PBH/data/mdot_table.h5'
    f = h5py.File(TabFile, 'r')
    mdot_table = f['mdot_table'][:]
    m_axis = f['m_axis'][:]
    z_axis = f['z_axis'][:]
    Light_BH_Table = {'mdot_table' : mdot_table,
                      'm_axis' : m_axis,
                      'z_axis' : z_axis}
    return Light_BH_Table

# Use as global variable to speed up calculations
Light_BH_Table = Get_BH_Interp_Data() 

def get_mdot(z = 10,
             m = 100):
    '''
    A fast interface for mdot,
    use interpolation for low masses to speed up
    '''
    m_axis = Light_BH_Table['m_axis']
    mdot_table = Light_BH_Table['m_amdot_tablexis']
    
    zp_axis = Light_BH_Table['z_axis'] + 1
    zp = z + 1
    
    if Within_Range(m, m_axis) and Within_Range(z, z_axis):
        # Use interpolation
        r = Interp_2D(
            Tab = mdot_table, 
            x_axis = m_axis, 
            y_axis = zp_axis,
            x_target = m,
            y_target = zp,
            Use_Log_X = True,
            Use_Log_Y = True,
            Use_Log_Z = True,
            )
        return r
    
    # If not in range then following lines will be executed
    # the inputs bellow must conform with interpolation table
    r = mdot(
        z = z,
        M_PBH = m,
        Omega_b = 0.0484,
        Use_halo = 1,
        Use_feedback = 0
    )
    return r
