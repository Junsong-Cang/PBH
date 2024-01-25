reload = 0
T21 = -100
datafile = 'data/6.npz'
ncpu = 12
nk = 50
Sigma = 0.5
m1 = 1
m2 = 1e3

LineWidth = 2
FontSize = 18

# Initialising
from src.Recombination import *
k1 = m2k(m2)
k2 = m2k(m1)
kc_vec = np.logspace(np.log10(k1), np.log10(k2), nk)

def fun(x):
    r = Find_PsMax_EDGES(
        T21 = T21,
        kc = x,
        SigmaPs = Sigma,
        LgF_Min = -30,
        zmax = 500,
        nz = 1500,
        nm = 500,
        DeltaC = 0.45,
        dmdm0_method = 0,
        Use_Halo_Boost = False,
        show_status = False,
        Use_Halo_Boost_Interp = True,
        Use_Halo_Boost_Extrap = True,
        Use_mdot_Interp = True,
        Use_LowM_Extrap_mdot = True,
        Use_HighM_Extrap_mdot = True,
        Fix_mdot = False,
        Fixed_mdot = 10,
        Show_Extrap_MSG = False,
        Ignore_Evo = False)
    SaySomething('6.log')
    return r

if reload:
    t1 = TimeNow()
    r = Parallel(n_jobs = ncpu)(delayed(fun)(x) for x in kc_vec)
    Timer(t1)
    np.savez(datafile, r = r)

p = np.load(datafile)['r']
bad_idx = []
for idx in np.arange(0, nk):
    if np.isnan(p[idx]):
        bad_idx.append(idx)

k = np.delete(kc_vec, bad_idx)
p = np.delete(p, bad_idx)

# Load OmB
r0 = np.load('data/3.npz')
k0 = r0['kc_vec']
p0 = r0['r']
print(p0)

bad_idx = []
for idx in np.arange(0, len(k0)):
    if np.isnan(p0[idx]):
        bad_idx.append(idx)

k0 = np.delete(k0, bad_idx)
p0 = np.delete(p0, bad_idx)

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

plt.plot(k0, p0, 'k', linewidth=LineWidth, label = '$\delta \Omega_{\mathrm{b}}/\Omega_{\mathrm{b}} < 30\%$')
plt.plot(k, p, 'r', linewidth=LineWidth, label = '$T_{21} < -100{\mathrm{mK}}$')

plt.xscale('log')
# plt.yscale('log')

plt.xlabel('$k\ [{\mathrm{Mpc^{-1}}}]$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$P_{\mathrm{s}}$',fontsize=FontSize,fontname='Times New Roman')
plt.legend(fontsize=FontSize,loc = 'upper left')
plt.ylim([0, 0.1])

plt.title('$\sigma = 0.5$',fontsize=FontSize)

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png',bbox_inches='tight',dpi=1000)
