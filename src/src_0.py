import numpy as np
from PyLab import *

def Vrel(z):
    '''
    Relative velosity in m/s
    See Eq.(3) of 2108.11130
    '''
    r = 3E4 * min(1, (1+z)/1e3)
    return r

def Cs(z):
    '''
    Sound speed in m/s
    See Eq.(2) of 2108.11130
    '''
    # Bolzmann constant
    kB = 1.38064852E-23
    # Proton mass in kg
    mp = 1.67262158E-27
    y = 5/3
    xe, Tk = LCDM_HyRec(z = z)
    r = np.sqrt(y * (1+xe) * Tk * kB/mp)
    return r

def Veff(z = 10.0, model = 1):
    '''
    ----inputs----
    model : veff model
    '''
    vr = Vrel(z)
    cs = Cs(z)
    if model == 0:
        r = np.sqrt( vr**2 + cs**2 )
    else:
        # Zhang's model
        if cs > vr:
            r = cs
        else:
            r = np.sqrt(cs*vr)
    return r

def Beta_naked(m, z, veff_model = 1):
    '''
    Viscosity for naked pbh
    See Eq.(3) of 2003.12589
    ----inputs----
    m : bh mass in msun
    '''
    veff = Veff(z = z, model = veff_model)
    xe, tmp = LCDM_HyRec(z = z)
    b1 = (m/1E4) * pow((1+z)/1e3, 1.5) * pow(5740/veff, 3)
    b2 = 0.257 + 1.45 * (xe/0.01) * pow((1+z)/1e3, 2.5)
    r = b1*b2
    return r

def Lambda_naked(m, z, veff_model = 1):
    '''
    ----inputs----
    m : pbh mass in msun
    '''
    b = Beta_naked(m = m, z = z, veff_model = veff_model)
    xcr = ( (1+b)**0.5 - 1 ) / b
    r = xcr**2 * np.exp(4.5/(3 + b**0.75))
    return r

def mdot_naked(m, z, veff_model = 1):
    '''
    Dimensionless accretion rate
    See Eq.(6) of 2108.11130
    ----inputs----
    m : pbh mass in msun
    '''
    l = Lambda_naked(m = m, z = z, veff_model = veff_model)
    veff = Veff(z = z, model = veff_model)
    r = 0.4 * l * ((1+z)/1e3)**3 * m * (1e3/veff)**3
    return r

def K_factor(m, z, veff_model = 1):
    '''
    K factor in Eq.(7) of 2108.11130
    '''
    mh = 3 * m * 1000 / (1+z)
    veff = Veff(z = z, model = veff_model)
    r = 0.22 * ((1+z)/1000) * mh**(2/3) * (1000/veff)**2
    return r

def Lambda_clothed(m, z, veff_model = 1, xcr_model = 0, xcr_model2 = 0):
    k = K_factor(m = m, z = z, veff_model = veff_model)
    b0 = Beta_naked(m = m, z = z, veff_model = veff_model)
    p = 0.75
    p1p = p/(1-p)
    bh = k**p1p * b0
    y = (1 + 10*bh)**0.1 * np.exp(2-k) * (k/2)**2
    
    # xcr, this is where things might go wrong
    if xcr_model == 0:
        b = b0
    else:
        b = bh
    x0 = ((1+b)**0.5 - 1)/b
    if xcr_model2 == 0:
        xcr = x0
    else:
        xcr = (k/2)**p1p * x0
    lbh = xcr**2 * np.exp(4.5/(3+bh**0.75))
    r = y**p1p * lbh
    return r

def mdot_clothed(m, z, veff_model = 1, xcr_model = 0, xcr_model2 = 0, Forbid_super_Eddington = 0):
    k = K_factor(m=m, z=z, veff_model=veff_model)
    if k>=2:
        mh = 3*m*1e3/(1+z)
        r = mdot_naked(m=m,z=z,veff_model=veff_model)
    else:
        L = Lambda_clothed(m=m,z=z,veff_model=veff_model,xcr_model=xcr_model, xcr_model2=xcr_model2)
        ve = Veff(z=z, model=veff_model)
        r = 3 * L * (1+z)/1e3 * m *(1e3/ve)**3
    if Forbid_super_Eddington:
        # need to account for radiation efficiency so maximum is not 1
        r = min(10, r)
    return r

def Mdot_Eddington(m):
    '''
    Eddington accretion rate, in kg/s
    '''
    r = 1.44E14 * m
    return r

def dmdlna(m, lna, veff_model = 1, xcr_model = 0, xcr_model2 = 0, Forbid_super_Eddington = 0):
    '''
    dm/dlna in msun unit
    '''
    msun = 2e30
    z = np.exp(-lna) - 1
    mdot_unitless = mdot_clothed(m =m, z=z, veff_model=veff_model, xcr_model=xcr_model,xcr_model2 = xcr_model2,Forbid_super_Eddington=Forbid_super_Eddington)
    dMdt = mdot_unitless * Mdot_Eddington(m) / msun
    r = dMdt/Hubble(z)
    return r

def Mbh_Evolution(m0 = 1e3,
                  zmin = 7,
                  zmax = 2e3,
                  nz = 10000,
                  veff_model = 1,
                  xcr_model = 1,
                  xcr_model2 = 1,
                  Forbid_super_Eddington = 0,
                  return_final_m = False
                  ):
    '''
    Evolve PBH mass
    '''
    x1 = np.log(1/(1+zmax))
    x2 = np.log(1/(1+zmin))
    y1 = m0
    f = lambda x, m: dmdlna(m,x,veff_model,xcr_model,xcr_model2,Forbid_super_Eddington)
    x, m = Integrate(f,y1,x1,x2,nx=nz)
    z = np.exp(-x) - 1
    if return_final_m:
        return m[-1]
    else:
        return z, m

'''
reload = 1
nm = 50
nz = 10000
zmin = 7
zmax = 2e3
m = np.logspace(0,3,nm)

f = np.empty((16,nm))

if reload:
    model = 0
    import time
    start = time.time()
    for veff_model in [0,1]:
        for xcr1 in [0,1]:
            for xcr2 in [0,1]:
                for edd in [0,1]:
                    for idx in np.arange(0,nm):
                        f[model,idx] = Mbh_Evolution(m[idx],zmin,zmax,nz,veff_model,xcr1,xcr2,edd,True)/m[idx]
                        # f[model,idx] = 1
                    model = model+1
                    print('model = ', model, ', time used :', time.time()-start)
    np.savez('tmp.npz',f=f)
else:
    f = np.load('tmp.npz')['f']


LineWidth = 2
FontSize = 18

import matplotlib.pyplot as plt
import os
plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
for idx in np.arange(0,16):
    plt.loglog(m,f[idx,:],color = 'k',linestyle='-',linewidth=LineWidth,label = str(idx))
# plt.legend(fontsize=FontSize,loc = 'lower left')
plt.xlabel('$m_0 [m_{\odot}]$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$m_f/m_0$',fontsize=FontSize,fontname='Times New Roman')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png',bbox_inches='tight',dpi=1000)
os.system('open /Users/cangtao/Desktop/tmp.png')
'''

'''
reload = 1
m1 = 1e0
m2 = 1e3
nz = 100
z = np.logspace(0,3,nz)
f1 = np.linspace(0,1,nz)
f2 = np.linspace(0,1,nz)
LineWidth = 2
FontSize = 18

for zid in np.arange(0,nz):    
    f1[zid] = mdot_clothed(m1,z[zid],1,1,1,0)
    f2[zid] = mdot_clothed(m2,z[zid],1,1,1,0)
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
plt.loglog(z,f1,color = 'k',linestyle='-',linewidth=LineWidth,label = 'm1')
plt.loglog(z,f2,color = 'r',linestyle='-',linewidth=LineWidth,label = 'm2')
plt.legend(fontsize=FontSize,loc = 'lower left')
plt.xlim([1, 1e3])
plt.ylim([1e-7,1e2])
plt.tight_layout()

plt.show()
'''



reload = 1
nm = 50
nz = 10000
zmin = 7
zmax = 2e3
m = np.logspace(0,3,nm)

f = np.linspace(0,1,nm)

for idx in np.arange(0,nm):
    f[idx] = Mbh_Evolution(m[idx],zmin,zmax,nz,1,1,1,1,True)/m[idx]
    print(idx/nm)

LineWidth = 2
FontSize = 18

import matplotlib.pyplot as plt
import os
plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
plt.loglog(m,f,color = 'k',linestyle='-',linewidth=LineWidth,label = 'test')
# plt.legend(fontsize=FontSize,loc = 'lower left')
plt.xlabel('$m_0 [m_{\odot}]$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$m_f/m_0$',fontsize=FontSize,fontname='Times New Roman')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png',bbox_inches='tight',dpi=1000)
os.system('open /Users/cangtao/Desktop/tmp.png')
