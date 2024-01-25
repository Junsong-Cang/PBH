
reload = 0
nm = 100
datafile = 'data/1.npz'
OmB_Lost_Ratio = 0.3
ncpu = 11

LineWidth = 2
FontSize = 18

from src.main import *

m = np.logspace(0, 5, nm)

def model(x, EoR, HaloBoost):

    r = Find_fbhmax_lite(
            m0 = x,
            Use_EoR = EoR,
            Use_Halo_Boost = HaloBoost,
            OmB_Lost_Ratio = OmB_Lost_Ratio,
            OmB0 = cosmo['OmB'],
            zmin = 8,
            zmax = 1000,
            nz = 2000,
            lf_min = -30,
            Precision = 1e-2,
            Use_Halo_Boost_Interp = True,
            Use_Halo_Boost_Extrap = True,
            Use_mdot_Interp = True,
            Use_LowM_Extrap_mdot = True,
            Use_HighM_Extrap_mdot = True,
            Show_Extrap_MSG = False)

    SaySomething('tmp.txt')
    
    return r

if reload:
    t1 = TimeNow()
    f = Parallel(n_jobs = ncpu)(delayed(model)(x, EoR = True, HaloBoost = False) for x in m)
    f1 = Parallel(n_jobs = ncpu)(delayed(model)(x, EoR = False, HaloBoost = False) for x in m)
    f2 = Parallel(n_jobs = ncpu)(delayed(model)(x, EoR = True, HaloBoost = True) for x in m)
    Timer(t1)
    np.savez(datafile, f = f, f1 = f1, f2 = f2, m = m)

r = np.load(datafile)
f = r['f']
f1 = r['f1']
f2 = r['f2']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

plt.loglog(m, f, 'k', linewidth=LineWidth, label = 'No Boost')
plt.loglog(m, f2, 'b', linewidth=LineWidth, label = 'Use Boost')

plt.xlabel('$m_0\ [m_{\odot}]$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$f_{\mathrm{bh}}$',fontsize=FontSize,fontname='Times New Roman')
plt.legend(fontsize=FontSize,loc = 'upper right')
plt.title('$\delta \Omega_{\mathrm{b}}<30 \%$',fontsize=FontSize)

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png',bbox_inches='tight',dpi=1000)
