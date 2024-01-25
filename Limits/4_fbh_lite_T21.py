from src.Recombination import *

reload = 0
T21 = -100
nm = 200
lm1 = -1
lm2 = 4
ncpu = 10
datafile = 'data/4.npz'
LineWidth = 2
FontSize = 18

m = np.logspace(lm1, lm2, nm)

def model(x, Ignore_Evo):
    r= Find_FbhMax_EDGES(
        T21 = T21,
        mc = x,
        sbh = -1,
        zmax = 1000,
        nz = 2000,
        nm = 1000,
        LgF_Min = -20,
        dmdm0_method = 0,
        Use_Halo_Boost = False,
        sbh_range = [-3, 5],
        show_status = False,
        Use_Halo_Boost_Interp = True,
        Use_Halo_Boost_Extrap = True,
        Use_mdot_Interp = True,
        Use_LowM_Extrap_mdot = True,
        Use_HighM_Extrap_mdot = True,
        Fix_mdot = False,
        Fixed_mdot = 10,
        Show_Extrap_MSG = False,
        Ignore_Evo = Ignore_Evo)
    SaySomething('4.log')
    return r

if reload:
    t1 = TimeNow()
    f = Parallel(n_jobs=ncpu)(delayed(model)(x, 0) for x in m)
    f1 = Parallel(n_jobs=ncpu)(delayed(model)(x, 1) for x in m)
    Timer(t1)
    np.savez(datafile, f1 = f1, f = f, m = m)

r = np.load(datafile)
f = r['f']
f1 = r['f1']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

plt.loglog(m, f, 'k', linewidth=LineWidth, label = 'Fiducial')
plt.loglog(m, f1, 'r', linewidth=LineWidth, label = 'No Mass Evo')

plt.xlabel('$m_0\ [m_{\odot}]$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$f_{\mathrm{bh}}$',fontsize=FontSize,fontname='Times New Roman')
plt.legend(fontsize=FontSize,loc = 'upper right')
plt.title('$T_{21} < -100 {\mathrm{mK}}$',fontsize=FontSize)

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png',bbox_inches='tight',dpi=1000)
