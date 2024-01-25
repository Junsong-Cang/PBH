reload = 0
datafile = 'data/9.npz'
ncpu = 12
nk = 40
ns = 40
m1 = 1e-1
m2 = 2e3
s1 = 0.2
s2 = 1.2
OmB_Lost_Ratio = 0.3

LineWidth = 2
FontSize = 18

# Initialising
from src.main import *
t = 318
t_tot = t * ns * nk/3600/12
print('Estimated time (hour): ', t_tot)

s_vec = np.linspace(s1, s2, ns)
k1 = m2k(m2)
k2 = m2k(m1)
kc_vec = np.logspace(np.log10(k1), np.log10(k2), nk)

N = ns*nk
idx_vec = np.arange(0, N)
params = np.zeros((N,2))

idx = 0
for kid in np.arange(0, nk):
    for sid in np.arange(0, ns):
        params[idx,0] = kc_vec[kid]
        params[idx,1] = s_vec[sid]
        # print(params[idx,:])
        idx = idx + 1
        
def fun(id):
    
    kc = params[id, 0]
    s = params[id, 1]
    
    r = Find_PsMax(
            kc = kc,
            Sigma = s,
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
    SaySomething('9.log')
    return r

if reload:
    t1 = TimeNow()
    r = Parallel(n_jobs = ncpu)(delayed(fun)(x) for x in idx_vec)
    Timer(t1)
    np.savez(datafile, r = r)
r = np.load(datafile)['r']
p = np.zeros((ns, nk))

idx = 0
for kid in np.arange(0, nk):
    for sid in np.arange(0, ns):
        p[sid,kid] = r[idx]
        idx = idx + 1

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

k, s = np.meshgrid(kc_vec, s_vec)

fig,ax = plt.subplots()
c=ax.pcolor(k, s, p, cmap='jet')#,norm = LogNorm(vmin=f.min(), vmax=f.max()))
plt.xlabel('$k\ [{\mathrm{Mpc^{-1}}}]$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$\sigma$',fontsize=FontSize,fontname='Times New Roman')
plt.xscale('log')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
clb = plt.colorbar(c)

plt.title('$P_{\mathrm{S}}$',fontsize=FontSize)

plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.png',dpi=200)
