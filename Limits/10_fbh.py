reload = 0
nm = 60
ns = 60
m1 = 1e0
m2 = 5e2
s1 = 0.1
s2 = 1.02
datafile = 'data/10.npz'
OmB_Lost_Ratio = 0.3
ncpu = 12

LineWidth = 2
FontSize = 18

from src.main import *
from matplotlib.colors import LogNorm

N = ns*nm
t = 80
T = t*N/3600/12
print('Estimated time (hour):', T)

lm1, lm2 = np.log10(m1), np.log10(m2)
m = np.logspace(lm1, lm2, nm)
s = np.linspace(s1, s2, ns)

idx_vec = np.arange(0, N)
params = np.zeros((N,2))

idx = 0
for sid in np.arange(0, ns):
    for mid in np.arange(0, nm):
        params[idx,0] = m[mid]
        params[idx,1] = s[sid]
        # print(params[idx,:])
        idx = idx + 1

def model(id):
    mc = params[id, 0]
    sbh = params[id, 1]
    
    r =  Find_fbhmax(
        mc = mc,
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
    
    SaySomething('10.txt')
    
    return r

if reload:
    t1 = TimeNow()
    f = Parallel(n_jobs = ncpu)(delayed(model)(x) for x in idx_vec)
    Timer(t1)
    np.savez(datafile, f = f)

r = np.load(datafile)['f']
f = np.zeros((ns, nm))

idx = 0
for sid in np.arange(0, ns):
    for mid in np.arange(0, nm):
        f[sid,mid] = r[idx]
        # print(params[idx,:])
        idx = idx + 1

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

m, s = np.meshgrid(m, s)

fig,ax = plt.subplots()
c=ax.pcolor(m, s, f,cmap='jet',norm = LogNorm(vmin=f.min(), vmax=f.max()))

plt.xlabel('$m_c\ [m_{\odot}]$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$\sigma$',fontsize=FontSize,fontname='Times New Roman')
plt.xscale('log')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
clb = plt.colorbar(c)

plt.title('$f_{\mathrm{bh}}$',fontsize=FontSize)

plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.png',dpi=200)
