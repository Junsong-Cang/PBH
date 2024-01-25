
reload = 0
nm = 50
sbh = 1
datafile = 'data/2.npz'
OmB_Lost_Ratio = 0.3
ncpu = 11

LineWidth = 2
FontSize = 18

from src.main import *

m = np.logspace(0, 4, nm)

def model(x):
    r =  Find_fbhmax(
        mc = x,
        sbh = sbh,
        OmB_Lost_Ratio = OmB_Lost_Ratio,
        Use_EoR = True,
        Use_Halo_Boost = False,
        sbh_range = [-3, 5],
        zmin = 8,
        zmax = 500,
        nz = 1500,
        nm = 500,
        OmB0 = cosmo['OmB'],
        lf_min = -30,
        Precision = 1e-2,
        show_status = False,
        Use_Halo_Boost_Interp = True,
        Use_Halo_Boost_Extrap = True,
        Use_mdot_Interp = True,
        Use_LowM_Extrap_mdot = True,
        Use_HighM_Extrap_mdot = True,
        Show_Extrap_MSG = False,
        Fix_mdot = False,
        Fixed_mdot = 10)
    
    SaySomething('2.txt')
    
    return r

if reload:
    t1 = TimeNow()
    f = Parallel(n_jobs = ncpu)(delayed(model)(x) for x in m)
    Timer(t1)
    np.savez(datafile, f = f)

r = np.load(datafile)
f = r['f']

r0 = np.load('data/1.npz')
m0 = r0['m']
f0 = r0['f']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

plt.loglog(m0, f0, 'k', linewidth=LineWidth, label = '$\sigma_{\mathrm{bh}} = 0$')
plt.loglog(m, f, 'r', linewidth=LineWidth, label = '$\sigma_{\mathrm{bh}} = 1$')

plt.xlabel('$m_0\ [m_{\odot}]$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$f_{\mathrm{bh}}$',fontsize=FontSize,fontname='Times New Roman')
plt.legend(fontsize=FontSize,loc = 'upper right')
plt.title('$\delta \Omega_{\mathrm{b}}<30 \%$',fontsize=FontSize)

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png',bbox_inches='tight',dpi=1000)
