reload = 0
datafile = 'data/3.npz'
ncpu = 12
nk = 50
Sigma = 0.5
m1 = 1
m2 = 1e3
OmB_Lost_Ratio = 0.3

LineWidth = 2
FontSize = 18

# Initialising
from src.main import *
k1 = m2k(m2)
k2 = m2k(m1)
kc_vec = np.logspace(np.log10(k1), np.log10(k2), nk)

def fun(x):
    r = Find_PsMax(
            kc = x,
            Sigma = Sigma,
            OmB_Lost_Ratio = OmB_Lost_Ratio,
            Use_EoR = True,
            DeltaC = 0.45,
            map_nx = 500,
            OmB0 = cosmo['OmB'],
            zmin = 8,
            zmax = 500,
            nz = 1500,
            nm = 500,
            Use_Halo_Boost = False,
            sbh_range = [-3, 5],
            show_status = False,
            Use_Halo_Boost_Interp = True,
            Use_Halo_Boost_Extrap = True,
            Use_mdot_Interp = True,
            Use_LowM_Extrap_mdot = True,
            Use_HighM_Extrap_mdot = True,
            Precision = 1e-2,
            Show_Extrap_MSG = False,
            Fix_mdot = 0,
            Fixed_mdot = 10)
    SaySomething('3.log')
    return r

if reload:
    t1 = TimeNow()
    r = Parallel(n_jobs = ncpu)(delayed(fun)(x) for x in kc_vec)
    Timer(t1)
    np.savez(datafile, r = r, kc_vec = kc_vec)

p = np.load(datafile)['r']

bad_idx = []
for idx in np.arange(0, nk):
    if np.isnan(p[idx]):
        bad_idx.append(idx)

k = np.delete(kc_vec, bad_idx)
p = np.delete(p, bad_idx)

print(len(k))

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

plt.plot(k, p, 'k', linewidth=LineWidth, label = '$\sigma_{\mathrm{bh}} = 0$')
plt.xscale('log')
# plt.yscale('log')

plt.xlabel('$k\ [{\mathrm{Mpc^{-1}}}]$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$P_{\mathrm{s}}$',fontsize=FontSize,fontname='Times New Roman')
# plt.legend(fontsize=FontSize,loc = 'upper right')
plt.ylim([0, 0.1])

plt.title('$\delta \Omega_{\mathrm{b}} / \Omega_{\mathrm{b}} < 30 \%,\ \sigma = 0.5$',fontsize=FontSize)

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png',bbox_inches='tight',dpi=1000)
