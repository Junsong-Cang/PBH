type = 1
reload = 0
Use_EoR = 1
nz = 2000

LineWidth = 2
FontSize = 15

from src.main import *
from matplotlib.colors import LogNorm

fbh0_template = '/Users/cangtao/cloud/GitHub/PBH/data/1.npz'
fbh0_data = np.load(fbh0_template)
m0_vec = fbh0_data['m']
nm = len(m0_vec)

if Use_EoR:
    fbh0_vec = fbh0_data['f']
else:
    fbh0_vec = fbh0_data['f1']

def rbh_pc(fbh0, m0, z):
    RhoC_cmv = 2.775e11 * 0.11933 # msun/Mpc^3
    nbh = fbh0 * RhoC_cmv/m0*(1+z)**3 # 1/Mpc^3
    vbh = 1/nbh
    # averaged radius occupied by bh, in pc
    rbh = (3*vbh/(4*np.pi))**(1/3)*1e6
    return rbh

def import_rate(z = 10, m = 10, Use_EoR = 1):
    
    # accretion rate
    mdot_edd = 1.44e14 * m # edd rate in kg/s
    m_dot = Get_mdot(z = z, m = m, Use_EoR = Use_EoR)
    dmdt = m_dot * mdot_edd * 0.9 # kg/s

    # import dmdt
    OmBh2 = 0.02242
    RhoB = 1.879e-26 * OmBh2 * (1+z)**3 # kg/m^3
    v = Effective_Speed(z = z, Use_EoR = Use_EoR)
    r = Bondi_Radius(z = z, M_PBH = m, Use_EoR = Use_EoR) * 3.086E16

    dmdt_import = RhoB*v*r**2
    
    ratio = dmdt_import / dmdt

    return ratio

def bondi_ratio(m0, f0_method):
    if f0_method == 0:
        f0 = 0.8
    else:
        f0 = np.interp(x = np.log(m0), xp = np.log(m0_vec), fp = fbh0_vec)
    Evo = Get_Evo_Lite(fbh0 = f0, m0 = m0, nz = nz, Use_EoR = Use_EoR)
    z_vec = Evo['z']
    m_vec = m0 * Evo['m_ratio']
    OmB_ratio = Evo['OmB_ratio']
    r_ratio = np.zeros(nz)
    import_ratio = np.zeros(nz)

    for idx in np.arange(0, nz):
        z_ = z_vec[idx]
        m_ = m_vec[idx]
        r_bondi = Bondi_Radius(z = z_, M_PBH = m_, Use_EoR = Use_EoR)
        rbh = rbh_pc(fbh0 = f0, m0 = m0, z = z_)
        r_ratio[idx] = r_bondi/rbh
        import_ratio[idx] = import_rate(z =z_, m = m_, Use_EoR = Use_EoR)
    
    result = {'z_vec' : z_vec, 'r_ratio' : r_ratio, 'OmB_ratio':OmB_ratio, 'import_ratio': import_ratio}
    
    return result

if reload:
    
    t1 = TimeNow()
    
    bondi_0 = np.zeros((nz, nm))
    bondi_1 = np.zeros((nz, nm))
    OmB0 = np.zeros((nz, nm))
    OmB1 = np.zeros((nz, nm))
    IR0 = np.zeros((nz, nm))
    IR1 = np.zeros((nz, nm))
    
    for mid in np.arange(0, nm):
        m_ = m0_vec[mid]
        br0 = bondi_ratio(m0 = m_, f0_method = 0)
        br1 = bondi_ratio(m0 = m_, f0_method = 1)
        z0 = br0['z_vec']
        z1 = br1['z_vec']
        
        for zid in np.arange(0, nz):
            bondi_0[zid, mid] = br0['r_ratio'][zid]
            bondi_1[zid, mid] = br1['r_ratio'][zid]
            OmB0[zid, mid] = 1-br0['OmB_ratio'][zid]
            OmB1[zid, mid] = 1-br1['OmB_ratio'][zid]
            IR0[zid,mid] = br0['import_ratio'][zid]
            IR1[zid,mid] = br1['import_ratio'][zid]
            
    z = br0['z_vec']
    Timer(t1)
    np.savez('tmp.npz', z = z, bondi_0 = bondi_0, bondi_1 = bondi_1, OmB0 = OmB0, OmB1 = OmB1, IR0 = IR0, IR1 = IR1)
    
rst = np.load('tmp.npz')
z = rst['z']
bondi_0 = rst['bondi_0']
bondi_1 = rst['bondi_1']
OmB0 = rst['OmB0']
OmB1 = rst['OmB1']
IR0 = rst['IR0']
IR1 = rst['IR1']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

m, z = np.meshgrid(m0_vec, z)
PlotData = bondi_1
PlotData = OmB1
PlotData = IR0

fig,ax = plt.subplots()

c=ax.pcolor(m, z, PlotData, cmap='jet',norm = LogNorm(vmin=PlotData.min(), vmax=PlotData.max()))
#c=ax.pcolor(m, z, PlotData, cmap='jet')

plt.xlabel('$m_0$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$z$',fontsize=FontSize,fontname='Times New Roman')
plt.xscale('log')
plt.yscale('log')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

# Colorbar label setting is painstaking
# in many cases it's perhaps easier to just use title+colorbar
# cbar = plt.colorbar(c)
# cbar.ax.tick_params(labelsize=FontSize/1.2)
# ax.text(2.5E11,2.6,'$f_{\\mathrm{bh}}$',rotation=0, size=FontSize)

# fig.colorbar(c, ax=ax,label='$\sigma$')

clb = plt.colorbar(c)
#clb.ax.set_title('$r_{\mathrm{eff}}/r_{\mathrm{occupied}}$')

contours = plt.contour(m, z, PlotData, levels = [0.1, 1], colors=['w', 'k'])
plt.title('$r_{\mathrm{eff}}/r_{\mathrm{occupied}}$',fontsize=FontSize)
plt.title('$\delta \Omega_{\mathrm{b}}/\Omega_{\mathrm{b}}$',fontsize=FontSize)
plt.title('Import Ratio',fontsize=FontSize)

plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.png',dpi=200)

r = import_rate()
print(r)

