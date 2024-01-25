from PyLab import *

def cs_kernel(z):
    zdec = 130
    zp = 1+z
    r1 = 5700 * (zp/1000)**0.5
    b = 1.72
    r2 = 1 + ((1+zdec)/zp)**b
    r2 = r2**(-1/(2*b))
    r = r1*r2
    return r

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

def v_eff(z, veff_model):
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

def mdot_kernel(m, z, veff_model):
    zp = 1+z
    veff = v_eff(z, veff_model)
    mh = 3*m*1000/zp
    # mbh = mh
    mbh = m
    
    k = 0.22 * (zp/1000) * mh**(2/3) * (veff/1e3)**(-2)
    if k < 2:
        return np.nan
    xe, swap = LCDM_HyRec(z = z)
    cs = cs_kernel(z)
    b1 = (mbh/1e4) * (zp/1000)**1.5
    b2 = (veff/5740)**(-3)
    # b2 = (cs/5740)**(-3)
    
    b3 = 0.257 + 1.45 *(xe/0.01) * (zp/1000)**2.5
    b = b1 * b2 * b3

    xcr = ((1+b)**0.5 - 1)/b
    
    Lambda = np.exp(4.5/(3 + b**0.75)) * xcr**2
    
    r = 0.023 * Lambda * (zp/1000) * mbh * (veff/5740)**(-3)
    
    return r

