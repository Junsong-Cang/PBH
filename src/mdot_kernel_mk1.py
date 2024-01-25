# Get mdot following 0709.0524
from PyLab import *
import numpy as np

def Tk_ez(z):
    '''
    An analytic approximation for Tk assuming adiabatic cooling
    '''
    zp = 1+z
    b = 1.72
    p1 = 2730 * zp/1000
    OmBh2 = 0.02242
    zdec = 132 * (OmBh2/0.022)**0.4
    adec = 1/(1+zdec)
    
    a = 1/zp
    p2 = adec/((a**b + adec**b)**(1/b))
    r = p1 * p2
    return r

def get_veff(z, veff_model):
    zp = 1+z
    Tk = Tk_ez(z)
    cs = 5700 * (Tk/2730)**0.5
    DataPath = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/'
    
    if veff_model == 0:
        vr = min(1, zp/1000) * 30000
        veff = (vr**2 + cs**2)**0.5

    elif veff_model == 1:

        f = DataPath + '0709.0524.Fig2.vrel.txt'
        zp_axis, vrel_axis = Read_Curve(
            File = f,
            nx = 1000,
            model = 2,
            Convert_x = 1,
            Convert_y = 1
        )
        vrel_axis = vrel_axis * 1e3
        vr = np.interp(x = zp, xp = zp_axis, fp = vrel_axis)
        veff = (vr**2 + cs**2)**0.5

    elif veff_model == 2:
        f = DataPath + '0709.0524.Fig2.veff_A.txt'
        zp_axis, veff_axis = Read_Curve(
            File = f,
            nx = 1000,
            model = 2,
            Convert_x = 1,
            Convert_y = 1
        )
        veff_axis = veff_axis * 1e3
        veff = np.interp(x = zp, xp = zp_axis, fp = veff_axis)
    elif veff_model == 3:
        f = DataPath + '0709.0524.Fig2.veff_B.txt'
        zp_axis, veff_axis = Read_Curve(
            File = f,
            nx = 1000,
            model = 2,
            Convert_x = 1,
            Convert_y = 1
        )
        veff_axis = veff_axis * 1e3
        veff = np.interp(x = zp, xp = zp_axis, fp = veff_axis)

    elif veff_model == 4:
        f = DataPath + '0709.0524.Fig2.vrel.txt'
        zp_axis, vrel_axis = Read_Curve(
            File = f,
            nx = 1000,
            model = 2,
            Convert_x = 1,
            Convert_y = 1
        )
        vrel_axis = vrel_axis * 1e3
        vr = np.interp(x = zp, xp = zp_axis, fp = vrel_axis)

        mach = vr/cs

        if mach > 1:
            veff = cs * (16 * mach**3 /np.sqrt(2 * np.pi))**(1/6)
        else:
            veff = cs * (1+mach**2)**0.5
    
    elif veff_model == 5:
        f = DataPath + '0709.0524.Fig2.vrel.txt'
        zp_axis, vrel_axis = Read_Curve(
            File = f,
            nx = 1000,
            model = 2,
            Convert_x = 1,
            Convert_y = 1
        )
        vrel_axis = vrel_axis * 1e3
        vr = np.interp(x = zp, xp = zp_axis, fp = vrel_axis)
        
        mach = vr/cs

        if mach > 1:
            veff = cs * mach * (np.sqrt(2 / np.pi) * np.log(2 * mach/np.e))**(-1/3)
        else:
            veff = cs * (1+mach**2)**0.5

    return veff

def mdot_naked(m, z, xe, veff_model):

    zp = 1+z
    Tk = Tk_ez(z)
    cs = 5700 * (Tk/2730)**0.5
    veff = get_veff(z,veff_model)
    

    b1 = m/1e4 *(zp/1000)**1.5 * (cs/5740)**(-3)
    # b1 = m/1e4 *(zp/1000)**1.5 * (veff/5740)**(-3)
    
    b2 = 0.257 + 1.45 * (xe/0.01)*(zp/1000)**2.5
    b = b1 * b2

    x = ((1+b)**0.5 -1)/b
    Lambda = x**2 * np.exp(4.5/(3 + b**0.75))
    r = 1.8e-3 * Lambda * (zp/1000)**3 * m * (veff/5740)**(-3)
    
    return r

def beta0(m, z, xe, Tk):
    zp = 1+z
    cs = 5700 * (Tk/2730)**0.5
    b1 = m/1e4
    b2 = (zp/1000)**1.5
    b3 = (5740/cs)**3
    b4 = 0.257 + 1.45*(xe/0.01)*(zp/1000)**2.5
    r = b1*b2*b3*b4
    return r

def xcr0(beta):
    r = (-1 + (1+beta)**0.5)/beta
    return r

def mdot_kernel(m,z,xe,veff_model):
    zp = 1+z
    mh = 3*m*1000/zp
    Tk = Tk_ez(z)
    cs = 5700 * (Tk/2730)**0.5
    veff = get_veff(z, veff_model)
    
    x = 0.22 * mh**(2/3) * (zp/1000) * (1000/cs)**2
    # x = 0.22 * mh**(2/3) * (zp/1000) * (1000/veff)**2

    if x < 2:
        mbh = m
        p = 3-2.25
        p2 = p/(1-p)

        b0 = beta0(m=mbh, z=z, xe=xe, Tk=Tk)
        b = x**p2 * b0

        x0 = xcr0(b)
        xh = x0*(x/2)**p2

        L0 = xh**2 * np.exp(4.5/(3+b**0.75))
        y = (1+10*b)**0.1 * np.exp(2-x)*(x/2)**2
        Lambda = L0 * y**p2
        Lambda = np.nan
    else:
        mbh = mh
        b = beta0(m=mbh, z=z, xe=xe, Tk=Tk)
        xcr = xcr0(b)
        Lambda = xcr**2 * np.exp(4.5/(3+b**0.75))
    
    mdot = 0.016 * Lambda * (zp/1000) * mbh * (veff/5.74e3)**(-3)
    
    return mdot
