from src.main import *

LineWidth = 2
FontSize = 15
m1 = 5
m2 = 10
m3 = 100
reload = 1
m = np.logspace(0, 4, 200)
nz = 2000
EoR = 1

#----
from matplotlib.colors import LogNorm

nm = len(m)
if reload:
    t1 = TimeNow()
    mdot_vec = np.zeros((nz, nm))
    m_vec = np.zeros((nz, nm))
    for mid in np.arange(0,nm):
        r = Get_Evo_Lite(fbh0 = 1e-8, m0 = m[mid], nz = nz, Use_EoR = EoR)
        m_ratio = r['m_ratio']
        m_vec[:, mid] = m[mid] * m_ratio
        z = r['z']
        for zid in np.arange(0, nz):
            mdot_vec[zid,mid] = Get_mdot(z = z[zid], m = m_vec[zid,mid], Use_EoR = EoR, Use_Edd_Limit = 1)
    Timer(t1)
    np.savez('tmp.npz', z = z, m_vec = m_vec, mdot_vec = mdot_vec)

r = np.load('tmp.npz')
z = r['z']
m_vec = r['m_vec']
mdot_vec = r['mdot_vec']

mid1 = Find_Index(x = m1, x_axis = m)
mid2 = Find_Index(x = m2, x_axis = m)
mid3 = Find_Index(x = m3, x_axis = m)
m_vec_1 = m_vec[:, mid1]
m_vec_2 = m_vec[:, mid2]
m_vec_3 = m_vec[:, mid3]

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

m_axis = m_vec[-1,:]
m_axis, z = np.meshgrid(m_axis, z)

fig,ax = plt.subplots()
c=ax.pcolor(m_axis, z, mdot_vec, cmap='jet', norm = LogNorm(vmin=mdot_vec.min(), vmax=mdot_vec.max()))
plt.plot(m_vec_1, z, '--k', linewidth=LineWidth)
plt.plot(m_vec_2, z, '--k', linewidth=LineWidth)
plt.plot(m_vec_3, z, '--k', linewidth=LineWidth)

plt.xlabel('$m_{\mathrm{final}}\ [m_\odot]$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$z$',fontsize=FontSize,fontname='Times New Roman')
plt.xscale('log')
plt.yscale('log')
plt.ylim([8, 300])

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
clb = plt.colorbar(c)

contours = plt.contour(m_axis, z, mdot_vec, levels=[1, 9], colors=['w','k'])
plt.title('$\dot{m}$',fontsize=FontSize)

plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.png',dpi=200)
