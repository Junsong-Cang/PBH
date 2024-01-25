type = 1
reload = 0
FileName = '/Users/cangtao/Desktop/21cmFAST-data/EOS_2021.h5'
redshift = 0.1
max_redshift = 35
LC_Quantities = ('brightness_temp','Ts_box','xH_box','Tk_box','Trad_box')
GLB_Quantities = ('brightness_temp','Ts_box','xH_box','dNrec_box','z_re_box','Gamma12_box','J_21_LW_box','density','Trad_box','Tk_box','Fcoll')

LineWidth = 2
FontSize = 15

import py21cmfast as p21c
import time, os
from PyLab import *
from p21c_tools import *

user_params = p21c.UserParams(
  HII_DIM = 100,
  N_THREADS = 1,
  USE_RELATIVE_VELOCITIES = True,
  USE_INTERPOLATION_TABLES = True,
  FAST_FCOLL_TABLES = True,
  )

astro_params = p21c.AstroParams(
  F_STAR10 = -1.25,
  F_STAR7_MINI = -2.5,
  L_X = 40.5,
  L_X_MINI = 40.5,
  NU_X_THRESH = 500.0,
  t_STAR = 0.5,
  ALPHA_STAR = 0.5,
  ALPHA_STAR_MINI = 0,
  A_LW = 2.0,
  BETA_LW = 0.6,
  A_VCB = 1.0,
  BETA_VCB = 1.8,
  F_ESC10 = -1.35,
  F_ESC7_MINI = -1.35,
  ALPHA_ESC = -0.3,
  )

flag_options = p21c.FlagOptions(
  USE_MINI_HALOS = True,
  USE_MASS_DEPENDENT_ZETA = True,
  INHOMO_RECO = True,
  USE_TS_FLUCT = True,
)

# ---- Initialise ----
if reload:
  t1 = TimeNow()
  lc = p21c.run_lightcone(
    redshift=redshift, 
    max_redshift=max_redshift,
    astro_params=astro_params, 
    flag_options=flag_options,
    user_params = user_params,
    lightcone_quantities=LC_Quantities,
    global_quantities=GLB_Quantities
    )
  Timer(t1)
  try:
    os.remove(FileName)
  except:
    pass
  lc.save(FileName)

lc = p21c.LightCone.read(FileName)
z = lc.node_redshifts
tk = lc.global_Tk
xH = lc.global_xH
tb = lc.global_brightness_temp
print(np.max(z))

'''
# Compare my results with EOS21 paper, looks very good to me
xH_File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/2110.13919.Fig9.Black.txt'
tb_File = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/2110.13919.Fig11.Black.txt'

z2, xH2 = Read_Curve(xH_File)
z2, tb2 = Read_Curve(tb_File)

x = [tb, xH, tk]
yscale = ['linear', 'linear', 'log']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

fig, axs = plt.subplots(1, 3, sharex = False, sharey = False)
fig.set_size_inches(10, 3)

axs[0].plot(z, tb, 'k', linewidth = LineWidth, label='Me')
axs[0].plot(z2, tb2, 'r', linewidth = LineWidth, label='EOS21')
axs[0].legend(fontsize=FontSize/1.2, loc = 'lower right')
axs[0].set_title('$T_{21}$',fontsize=FontSize)
axs[0].set_xlabel('$z$',fontsize=FontSize,fontname='Times New Roman')
axs[0].tick_params(axis='both', which='both', labelsize = FontSize)

axs[1].plot(z, xH, 'k', linewidth = LineWidth, label='my codes')
axs[1].plot(z2, xH2, 'r', linewidth = LineWidth, label='EOS Paper')
axs[1].set_title('$x_{\mathrm{H}}$',fontsize=FontSize)
axs[1].set_xlabel('$z$',fontsize=FontSize,fontname='Times New Roman')
axs[1].tick_params(axis='both', which='both', labelsize = FontSize)
axs[1].set_ylim(0, 1.2)

axs[2].semilogy(z, tk, 'k', linewidth = LineWidth, label='my codes')
axs[2].set_title('$T_{\mathrm{k}}$',fontsize=FontSize)
axs[2].set_xlabel('$z$',fontsize=FontSize,fontname='Times New Roman')
axs[2].tick_params(axis='both', which='both', labelsize = FontSize)

plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=200)
'''

# Now joind 21cmFAST with LCDM HyRec

# HyRec Results
LCDM = '/Users/cangtao/cloud/GitHub/Dark_CosmoMC/HyRec/tmp_LCDM.dat'
d = np.loadtxt(LCDM)
z0 = d[:,0][::-1]
x0 = d[:,1][::-1]
t0 = d[:,2][::-1]


# 21cmFAST Results
xe = Switch_xe_format(xe = 1-xH, format = 1)

print(x0.max())
print(xe.max())

# xe = 1-xH
z = z[::-1]
tk = tk[::-1]
xe = xe[::-1]

# Get new templates, 1 EoR
nz0 = len(z0)
x1 = np.linspace(0,1,nz0)
t1 = np.linspace(0,1,nz0)

for zid in np.arange(0, nz0):
  
  z_ = z0[zid]

  if z_ > np.max(z):
    x1[zid] = x0[zid]
    t1[zid] = t0[zid]
  else:
    # ensure a smooth transition for xe
    x1[zid] = max(np.interp(x = z_, xp = z, fp = xe), np.interp(x = z_, xp = z0, fp = x0))
    t1[zid] = np.interp(x = z_, xp = z, fp = tk)

# Get them all to a shorter z-axis

nz = 200
zp = np.logspace(0, 3.4, nz)
x0_ = np.interp(x = zp, xp = z0 + 1, fp = x0)
t0_ = np.interp(x = zp, xp = z0 + 1, fp = t0)
x1_ = np.interp(x = zp, xp = z0 + 1, fp = x1)
t1_ = np.interp(x = zp, xp = z0 + 1, fp = t1)

lzp = np.log10(zp)
lx0_ = np.log10(x0_)
lx1_ = np.log10(x1_)
lt0_ = np.log10(t0_)
lt1_ = np.log10(t1_)
d = np.zeros((5, nz))
d[0,:] = lzp[:]
d[1,:] = lx0_[:]
d[2,:] = lt0_[:]
d[3,:] = lx1_[:]
d[4,:] = lt1_[:]

file = '/Users/cangtao/Desktop/tmp.txt'
np.savetxt(fname = file, X = d, fmt = '%.4f', delimiter = ', ')

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

fig, axs = plt.subplots(1, 2, sharex = False, sharey = False)
fig.set_size_inches(10, 4)

axs[0].loglog(1+z, xe, 'k', linewidth = LineWidth, label='21cmFAST')
axs[0].loglog(1+z0, x0, 'r', linewidth = LineWidth, label='LCDM')
axs[0].loglog(1+z0, x1, '-b', linewidth = LineWidth, label='Interpolated')
axs[0].legend(fontsize=FontSize/1.2, loc = 'lower right')
axs[0].set_title('$x_{\mathrm{e}}$',fontsize=FontSize)
axs[0].set_xlabel('$1+z$',fontsize=FontSize,fontname='Times New Roman')
axs[0].tick_params(axis='both', which='both', labelsize = FontSize)

axs[1].loglog(1+z, tk, 'k', linewidth = LineWidth, label='my codes')
axs[1].loglog(1+z0, t0, 'r', linewidth = LineWidth, label='LCDM')
axs[1].loglog(1+z0, t1, '--b', linewidth = LineWidth, label='ITP')
axs[1].set_title('$T_{\mathrm{k}}$',fontsize=FontSize)
axs[1].set_xlabel('$1+z$',fontsize=FontSize,fontname='Times New Roman')
axs[1].tick_params(axis='both', which='both', labelsize = FontSize)
# axs[1].set_ylim(0, 1.2)

plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=500)
