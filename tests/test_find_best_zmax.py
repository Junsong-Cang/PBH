# Find a stable set of zmax and nz (for zmin = 8)
# Stat: best precision for zmax = 500 is 2%
# nz = 1000, 1500, 2000 : 36%, 24%, 18%

reload = 0
LineWidth = 2
FontSize = 18

from src.main import *
if reload:
    r0 = Get_Evo_Lite(
        zmin = 8, 
        zmax = 1000,
        nz = 100000,
        Fix_mdot = 1,
        Fixed_mdot = 10
        )
    z = r0['z']
    m = r0['m_ratio']
    np.savez('tmp.npz', z = z, m = m)

r = np.load('tmp.npz')
z0 = r['z']
m0 = r['m']

# trick : find zmax first
r1 = Get_Evo_Lite(
    zmin = 8,
    zmax = 500,
    nz = 2000,
    Fix_mdot = 1,
    Fixed_mdot = 10
    )
z1 = r1['z']
m1 = r1['m_ratio']
m0_ = np.interp(x = z1, xp = z0[::-1], fp = m0[::-1])

dif = np.zeros(len(z1))
for idx in np.arange(0, len(z1)):
    M1 = m1[idx]
    M0 = m0_[idx]

    dif[idx] = np.abs(M1 - M0)/min(M1, M0)

dif_max = dif.max()
dif_avg = np.sum(dif)/len(dif)
print('dif_max = ', dif_max)
print('dif_avg = ', dif_avg)

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

plt.loglog(z0, m0, 'k', linewidth=LineWidth, label = 'Best')
plt.loglog(z1, m1, '--b', linewidth=LineWidth, label = 'nz = 1000')
plt.loglog(z0, 1.05*z0/z0, 'r', linewidth=LineWidth, label = '5\%')

plt.xlabel('$z$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$m/m_0$',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.legend(fontsize=FontSize,loc = 'upper right')
plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)

plt.show()
