import numpy as np
import h5py, os, time, warnings
from PyLab import *
import matplotlib.pyplot as plt
# from scipy.optimize import root as Sci_Solve

# convert between f[Hz] and k[Mpc^-1]
k2f = lambda k: k * 1.546e-15
f2k = lambda f: f / 1.546e-15

def Get_DoF(k = np.logspace(0, 20, 100)):
    '''
    Interpolate to find DoF, using template in
    /Users/cangtao/cloud/Library/Matlab/AUX/Cosmology.h5
    '''
    LgK_axis = np.linspace(2, 11, 20)
    LgDoF_axis = np.array([0.5271, 0.5233, 0.6073, 0.9083, 1.0140, 1.0228, 1.0246, 1.0387,
                1.1577, 1.3348, 1.7353, 1.8594, 1.9026, 1.9113, 1.9384, 1.9909, 2.0167, 2.0272, 2.0284, 2.0284])
    LgK_Min = np.min(LgK_axis)
    LgK_Max = np.max(LgK_axis)
    try:
        nk = len(k)
    except:
        k = np.array([k])
        nk = 1
    r = []
    for lk in np.log10(k):
        if lk < LgK_Min:
            DoF = 3.36
        elif lk > LgK_Max:
            DoF = 106.75
        else:
            DoF = 10**np.interp(lk, LgK_axis, LgDoF_axis)
        r.append(DoF)
    r = np.array(r)
    return r

def k2m(k):
    '''
    Get PBH mass in msun
    Let's just use the 2107.08638 results
    '''
    Gamma = 0.2
    DoF = Get_DoF(k)
    r = 2.43e12 * (Gamma/0.2) * (106.75/DoF)**(1/6) / k**2
    return r

def m2k(m):
    k_axis = np.logspace(-2, 25, 10000)
    m_axis = k2m(k_axis)
    # np.interp needs x_axis to be increasing
    m_axis = m_axis[::-1]
    k_axis = k_axis[::-1]
    r = np.interp(m, m_axis, k_axis)
    return r

def Power_Spectra(A = 1.0, kc = 1.0, Sigma = 1.0, k = np.logspace(-3, 3, 100)):
    '''
    Nomrlised PS model in NanoGrav15 New physics Paper
    '''
    x = np.log(k)
    xc = np.log(kc)
    y1 = ((x - xc)/Sigma)**2
    y2 = np.exp(-y1/2)
    Norm = np.sqrt(2*np.pi*Sigma**2)
    r = A*y2/Norm
    return r

# Some pre-requisites

def Integrand(A, kc, Sigma, k, v, u):
    '''
    Follow my own paper
    '''
    
    OmR = 9.1e-5
    h = 0.6766
    DoF = Get_DoF(k)
    # Easier for k to be a scalar
    if not Is_Scalar(k):
        raise Exception('k MUST be a scalar')

    # 1st line
    L1 = 0.29*OmR*h**2*(106.75/DoF)**(1/3)
    
    # 2nd line
    L2 = (4 * v**2 - (1 - u**2 + v**2)**2)/(4 * u**2 * v**2)
    L2 = L2**2

    # 3rd line
    F1a = np.log(abs((3 - (u+v)**2) / (3 - (u-v)**2)))
    F1b = 4*u*v/(u**2 + v**2 - 3)
    F1 = (F1a - F1b)**2
    F2 = np.pi**2 * np.heaviside(u + v - np.sqrt(3), 0)
    F = F1 + F2

    Pv = Power_Spectra(A = A, kc = kc, Sigma = Sigma, k = k*v)
    Pu = Power_Spectra(A = A, kc = kc, Sigma = Sigma, k = k*u)

    L3 = ((u**2 + v**2 - 3)/(2*u*v))**4
    L3 = L3 * F * Pv * Pu

    r = L1 * L2 * L3
    
    return r

def Integrate_u(A = 1.0, kc = 1.0, Sigma = 1.0, k = 1.0, v = 1.0, nu = 1000):
    
    fu = lambda x: Integrand(A = A, kc = kc, Sigma = Sigma, k = k, v = v, u = x)
    
    # Do a cut or you will get NaN error from numerical integration
    u1 = max(abs(1-v), 1e-4)
    u2 = 1+v
    
    u = np.linspace(u1, u2, nu)
    y = Function_Array(fu, u)
    r = np.trapz(y, u)
    return r

def Integrate_u_vec(A = 1.0, kc = 1.0, Sigma = 1.0, k = 1.0, v = np.linspace(0, 10, 20), nu = 1000):
    
    f = lambda x: Integrate_u(A = A, kc = kc, Sigma = Sigma, k = k, v = x, nu = nu)
    r = Function_Array(f, v)
    return r
    

def dOmGW0h2_dlnk_slow(A = 1.0, kc = 1.0, Sigma = 1.0, k = 1.0, nu = 1000, nv = 100):
    
    f = lambda v: Integrate_u_vec(A = A, kc = kc, Sigma = Sigma, k = k, v = v, nu = nu)
    v, fv = Map(
        F = f,
        Start = 1,
        Width = 0.6,
        MinX = -np.inf,
        MaxX = np.inf,
        nx = nv,
        Precision = 1e-2,
        Max_Iteration = 100,
        Use_Booster = 0
    )
    r = np.trapz(fv, v)
    return r

def dOmGW0h2_dlnk_slow_vec(A = 1.0, kc = 1.0, Sigma = 1.0, k = np.logspace(-3, -3, 50), nu = 1000, nv = 100, ncpu = 1):
    f = lambda x: dOmGW0h2_dlnk_slow(A = A, kc = kc, Sigma = Sigma, k = x, nu = nu, nv = nv)
    if ncpu == 1:
        r = Function_Array(f, k)
    else:
        r = Parallel(n_jobs=ncpu)(delayed(f)(x) for x in k)
    return r

def Convert_to_Primordial(k, GW0, Abort = True):
    '''
    Get primordial dOmegaGWh2/dlnk
    ----inputs----
    k : k
    GW0 : Current dOmegaGWh2/dlnk
    '''
    if Abort:
        # do nothing
        r = GW0
        return r

    OmR = 9.1e-5
    DoF = Get_DoF(k)
    Prefix = 0.38*OmR*(106.75/DoF)**(1/3)
    r = GW0 / Prefix
    return r

'''
# All above are the semi-analytic solutions, now try interpolation
def Get_Interp_Data():
    # This table contains dOmegaGW/dlnk for different Sigma and x (defined as k/kc)
    Small = 1e-200
    H5_File = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/GW_Tables.h5'
    f = h5py.File(H5_File, 'r')
    Sigma_axis = f['Sigma_axis'][:]
    x_axis = f['x_axis'][:]
    dGWdlnk_Tables = f['dGWdlnk_Tables'][:]
    # cutoff small number to avoid NaN error, can happen for inf - inf
    Size = np.shape(dGWdlnk_Tables)
    for id1 in np.arange(0, Size[0]):
        for id2 in np.arange(0, Size[1]):
            r = dGWdlnk_Tables[id1, id2]
            if r < Small:
                dGWdlnk_Tables[id1, id2] = Small
    Interp_Dataset = {'Sigma_axis':Sigma_axis, 'x_axis':x_axis, 'dGWdlnk_Tables':dGWdlnk_Tables}
    return Interp_Dataset

Interp_Dataset = Get_Interp_Data()
'''

def dOmGWh2_dlnk_delta(
        f = np.logspace(-4, 4, 100),
        fc = 1,
        k = np.logspace(-4, 4, 100),
        kc = 1.0,
        A = 1.0,
        Use_Freq_Domain = True,
        Use_Today_Value = True
        ):
    '''
    # Write your model here
    if fc < 1e200:
        raise Exception('Delta module not ready')
    '''
    h = 0.6766
    if Use_Freq_Domain:
        x = f/fc
        k_vec = f2k(f)
    else:
        x = k/kc
        k_vec = k
    
    x2 = x**2
    
    # The equation is very looooong
    p1 = 3*A**2/64
    p2 = ((4-x2)/4)**2
    p3 = x2
    p4 = (3*x2 - 2)**2
    L1 = p1*p2*p3*p4

    L2 = np.pi**2 * (3*x2 - 2)**2 * np.heaviside(2*np.sqrt(3) - 3*x, 0.5)

    L3a = (4 + (3*x2 - 2)*np.log(np.abs(1 - 4/(3*x2))))**2
    L3b = np.heaviside(2 - x, 0.5)

    r = L1*(L2 + L3a)*L3b*h**2

    if Use_Today_Value:
        OmR = 9.1e-5
        DoF = Get_DoF(k_vec)
        Prefix = 0.38*OmR*(106.75/DoF)**(1/3)
        r = Prefix * r

    return r

def dOmGWh2_dlnk(
        f = np.logspace(-4, 4, 100),
        fc = 1.0,
        k = np.logspace(-4, 4, 100),
        kc = 1.0,
        A = 1.0,
        Sigma = 0.5,
        Use_Freq_Domain = True,
        Use_Today_Value = True,
        Use_Delta_Approx = True,
        Use_Fast_Mode = True,
        nu = 1000,
        nv = 100,
        ncpu = 1
        ):
    '''
    Get dOmeGW0h2/dlnk
    ----inputs----
    f : frequency in Hz
    fc : central frequency in Hz
    k : k in Mpc^-1
    kc : kc (central) in Mpc^-1
    A : Amplitude
    Sigma : Width
    Use_Freq_Domain : Use frequency or k, unselected one will be ignored
    Use_Today_Value : Return current or primordial OmegaGW distribution
    '''

    '''
    Question:
    1. What if no solution found in interpolation table
    2. What if Sigma is too small: Use Delta
    '''

    if Use_Freq_Domain:
        x = f/fc
        kc_real = f2k(fc)
        k_vec = f2k(f)
    else:
        x = k/kc
        kc_real = kc
        k_vec = k

    if not Use_Fast_Mode:
            r = dOmGW0h2_dlnk_slow_vec(A = A, kc = kc_real, Sigma = Sigma, k = k_vec, nu = nu, nv = nv, ncpu = ncpu)
            r = Convert_to_Primordial(k = k_vec, GW0 = r, Abort = Use_Today_Value)
            return r

    # Prepare Interpolation data
    x_axis = Interp_Dataset['x_axis']
    Sigma_axis = Interp_Dataset['Sigma_axis']
    
    # Interpolaiton is done in log scale
    lx_axis = np.log10(x_axis)
    ls_axis = np.log10(Sigma_axis)
    Tab = np.log10(Interp_Dataset['dGWdlnk_Tables'])
    
    Sigma_Min = np.min(Sigma_axis)
    Sigma_Max = np.max(Sigma_axis)

    # Use Delta result if Sigma too small
    if Sigma < Sigma_Min:
        if Use_Delta_Approx:
            warnings.warn('Sigma too small for interpolation, using delta approximation')
            r = dOmGWh2_dlnk_delta(f = f, fc = fc, k = k, kc = kc, A = A, Use_Freq_Domain = Use_Freq_Domain, Use_Today_Value = Use_Today_Value)
            return r
        else:
            warnings.warn('Sigma too small for interpolation, numerical integration may fail')
            r = dOmGW0h2_dlnk_slow_vec(A = A, kc = kc_real, Sigma = Sigma, k = k_vec, nu = nu, nv = nv)
            r = Convert_to_Primordial(k = k_vec, GW0 = r, Abort = Use_Today_Value)
            return r
    elif Sigma > Sigma_Max:
            warnings.warn('Sigma too large for interpolation, numerical integration may fail')
            r = dOmGW0h2_dlnk_slow_vec(A = A, kc = kc_real, Sigma = Sigma, k = k_vec, nu = nu, nv = nv)
            r = Convert_to_Primordial(k = k_vec, GW0 = r, Abort = Use_Today_Value)
            return r
    
    # x can still have range overflow
    id1 = Find_Index(x = Sigma, x_axis = Sigma_axis)
    id2 = id1 + 1
    v = np.log10(Sigma)
    v1 = ls_axis[id1]
    v2 = ls_axis[id2]
    
    F1 = Tab[:,id1]
    F2 = Tab[:,id2]
    Fx_axis = (F2 - F1)*(v - v1)/(v2 - v1) + F1

    lx = np.log10(x)
    Fx0 = np.interp(x = lx, xp = lx_axis, fp = Fx_axis, left = np.nan, right = np.nan)
    Fx1 = 10**Fx0
    Fx2 = A**2 * Fx1 # re-scaled primordial value

    # Convert to doday
    OmR = 9.1e-5
    h = 0.6766
    DoF = Get_DoF(k_vec)
    Prefix = 0.38*OmR*h**2*(106.75/DoF)**(1/3)
    Fx = Fx2 * Prefix

    # Find overflow
    for idx in np.arange(0, len(Fx)):
        if np.isnan(Fx[idx]):
            warnings.warn('Found interpolation overflow, using slow mode')
            Fx[idx] = dOmGW0h2_dlnk_slow(A = A, kc = kc_real, Sigma = Sigma, k = k_vec[idx], nu = nu, nv = nv)
    
    r = Convert_to_Primordial(k = k_vec, GW0 = Fx, Abort = Use_Today_Value)
    
    return r

def Get_dNeff(
        fc = 1.0,
        kc = 1.0,
        A = 1.0,
        Sigma = 0.5,
        Use_Freq_Domain = False,
        nu = 1000,
        nv = 100,
        nx = 10000
        ):
    lxmin = np.log10(np.min(Interp_Dataset['x_axis']))
    lxmax = np.log10(np.max(Interp_Dataset['x_axis']))
    x_vec = np.logspace(lxmin, lxmax, nx)
    if Use_Freq_Domain:
        Kc = f2k(fc)
    else:
        Kc = kc
    k_vec = x_vec*Kc
    GWh2_vec = dOmGWh2_dlnk(
        k = k_vec,
        kc = Kc,
        A = A,
        Sigma = Sigma,
        Use_Freq_Domain = False,
        Use_Today_Value = True,
        Use_Delta_Approx = False,
        Use_Fast_Mode = True,
        nu = nu,
        nv = nv
        )
    h2 = 0.6766**2
    GW_vec = GWh2_vec/h2
    lnk = np.log(k_vec)
    GW = np.trapz(GW_vec, lnk)
    r = 8.3e4 * GW
    return r

def Variance_Kernel(
        A = 1,
        Sigma = 1,
        x = 1,
        y = 1
        ):
    '''
    Variance Kernel for PBH
    Further details can be found in Eq.5 of my SIGW paper.
    x = k/kc
    y = k'/k
    '''
    V_model = 0
    Prefix = 16*A/(81*Sigma*np.sqrt(2 * np.pi))
    f1 = (np.log(x*y))**2
    f2 = - y**2 - f1/(2*Sigma**2)
    f3 = y**3 * np.exp(f2)
    if V_model == 1:
        # Model in 2306.16219
        y3 = y/np.sqrt(3) 
        T = 3*(np.sin(y3) - y3*np.cos(y3))/y3**3
    else:
        T = 1
    r = Prefix*f3*T**2
    return r

def PBH_Variance(
        A = 1,
        Sigma = 1,
        x = np.logspace(-2, 2, 100),
        method = 0
        ):
    '''
    Speed : 2300 per second
    '''
    Is_Scalar(x, 1)
    nx = len(x)
    r = np.linspace(0, 1.0 ,nx)

    for idx in np.arange(0, nx):
        x_ = x[idx]
        f = lambda y: Variance_Kernel(A = A, Sigma = Sigma, x = x_, y = y)
        if method == 0:
            y_vec, f_vec = Map(
                F = f,
                Start = 0,
                Width = 0.3,
                MinX = -np.inf,
                MaxX = np.inf,
                nx = 500,
                Precision = 1e-3,
                Max_Iteration = 1200,
            )
        else:
            # This should give enough precision, if not then f**k it
            y_vec = np.linspace(-10, 10, 1000000)
            f_vec = f(y_vec)
        r[idx] = np.trapz(f_vec, y_vec)

    return r

def Phi_Kernel(
        A = 0.2,
        Sigma = 1.5,
        kc = 1e8,
        m = np.logspace(-8, 2, 1000),
        DeltaC = 0.45
        ):
    k = m2k(m)
    x = k/kc
    S2 = PBH_Variance(A = A, Sigma = Sigma, x = x)
    b1 = np.exp(-DeltaC**2 / (2*S2))
    b2 = np.sqrt(2*S2/(np.pi * DeltaC**2))
    Beta = b1*b2

    # Now get Phi
    Gamma = 0.2
    DoF = Get_DoF(k)
    r = 0.28 * (Beta*1e8) * (Gamma/0.2)**1.5 * (106.75/DoF)**0.25 / m**0.5

    return r

def LN_Profile_PBH(
        fbh = 1,
        Width = 1,
        mc = 1e4,
        m = np.logspace(-4, 10, 100)
    ):
    Is_Scalar(m,1)
    Prefix = fbh/(Width*np.sqrt(2*np.pi))
    p1 = (np.log(m/mc))**2
    p2 = - p1/(2*Width**2)
    p3 = np.exp(p2)
    r = Prefix*p3
    return r
    
def Phi_Profile(
        A = 0.2,
        Sigma = 1.5,
        kc = 1e8,
        DeltaC = 0.45,
        map_nx = 500
        ):
    Small = 1e-100
    Width_Min = 1e-2
    Width_Max = 20
    
    f = lambda m: Phi_Kernel(A = A, Sigma = Sigma, kc = kc, DeltaC = DeltaC, m = m)
    Peak = np.log10(k2m(kc))
    m_vec, Phi_vec = Map(
        F = f,
        Start = Peak,
        Width = 1,
        MinX = -np.inf,
        MaxX = np.inf,
        nx = map_nx,
        Precision = 1e-3,
        Max_Iteration = 1000,
        )
    
    lnm = np.log(m_vec)
    fbh = np.trapz(Phi_vec, lnm) #--------Need this
    if fbh < Small:
        Width = Width_Max #----Need this
        mc = 10**Peak #----Need this
        Best_Fit = LN_Profile_PBH(fbh = fbh, Width = Width, mc = mc, m = m_vec)
        # Ok let's raise errors deliberately
        r = {'fbh':fbh, 'mc' : np.nan, 'Width': np.nan, 'm_vec' : np.nan, 'Phi_vec' : np.nan, 'Best_Fit' : np.nan}
        return r
    
    Mc_idx = np.argmax(Phi_vec)
    if not Is_Scalar(Mc_idx):
        Mc_idx = Mc_idx[0]
    mc = m_vec[Mc_idx]
    
    # Now find Width
    dw = 0.01
    Phi_Theory = lambda w: LN_Profile_PBH(fbh = fbh, Width = w, mc = mc, m = m_vec)
    Chi2 = lambda w: np.sum((Phi_Theory(w) - Phi_vec)**2)
    dC2dw = lambda w: (Chi2(w+dw) - Chi2(w))/dw
    Width = Solve(
        F = dC2dw,
        Xmin = Width_Min,
        Xmax = 50,
        Precision = 1e-4
        )
    Best_Fit = LN_Profile_PBH(fbh = fbh, Width = Width, mc = mc, m = m_vec)
    r = {'fbh':fbh, 'mc':mc, 'Width': Width, 'm_vec':m_vec, 'Phi_vec':Phi_vec, 'Best_Fit':Best_Fit}
    return r
