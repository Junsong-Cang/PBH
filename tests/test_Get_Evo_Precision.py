# zmax, nz = 500, 1500 should be enough
# Stats : 
# [zmax, nz, nm] : Error@z10, Error@z8
# 500, 1500, 500 : 6%, None 
# 500, 1000, 500 : 12%, 20%
# 500, 1000, 200 : 12%, 20%
# Stangely nm = 200 shows no decrease in precision

reload = 1
datafile = 'tmp.npz'
LineWidth = 2
FontSize = 15

from src.main import *
from matplotlib.colors import LogNorm

if reload:
    t1 = TimeNow()
    r = Get_Evo(
        fbh0 = 1e-6,
        mc = 1e1,
        sbh = 2,
        nz = 4000,
        zmax = 1000,
        Fix_mdot = 0,
        Fixed_mdot = 10,
        nm = 1000,
        show_status = 1)
    Timer(t1)
    #np.savez(datafile, MF = r['MF'], m_vec = r['m_vec'], z_vec = r['z_vec'], 
    #         omb_ratio = r['omb_ratio'], m0_vec = r['m0_vec'], fbh_ratio = r['fbh_ratio'], m_ratio = r['m_ratio'])
    np.savez(datafile, z_vec = r['z_vec'], fbh_ratio = r['fbh_ratio'])

r0 = np.load(datafile)

z0 = r0['z_vec']
f0 = r0['fbh_ratio']

t1 = TimeNow()

r1 = Get_Evo(
    fbh0 = 1e-6,
    mc = 1e1,
    sbh = 2,
    nz = 1000,
    zmax = 500,
    Fix_mdot = 0,
    Fixed_mdot = 10,
    nm = 500,
    show_status = 1)

Timer(t1)
z1 = r1['z_vec']
f1 = r1['fbh_ratio']

F0 = np.interp(x = 10, xp = z0[::-1], fp = f0[::-1])
F1 = np.interp(x = 10, xp = z1[::-1], fp = f1[::-1])
dif = np.abs(F0-F1)/(min(F0, F1))
print(dif)

# Get plot
plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

plt.loglog(z0, f0, 'k', linewidth=LineWidth, label = 'HiRes')
plt.loglog(z1, f1, 'b', linewidth=LineWidth, label = 'LowRes')

plt.loglog(z0, 1.05*f0/f0, 'r', linewidth=LineWidth, label = '$5\%$')

plt.xlabel('$z$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('fbh ratio',fontsize=FontSize,fontname='Times New Roman')
plt.xscale('log')
plt.yscale('log')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.legend(fontsize=FontSize,loc = 'upper right')
plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=500)

plt.show()
