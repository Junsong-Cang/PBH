from src.Tau import *

reload = 0
m1 = 1e-3
m2 = 1e4
nm = 100
ncpu = 11
datafile = 'data/8.npz'

LineWidth = 2
FontSize = 18

m = np.logspace(np.log10(m1), np.log10(m2), nm)
def model(x, s):
    if s<0:
        r = Find_fbh_Tau(mc = x, sbh = s, zmax = 1000)
    else:
        r = Find_fbh_Tau(mc = x, sbh = s)
    SaySomething('8.log')
    return r

if reload:
    t1 = TimeNow()
    r0 = Parallel(n_jobs=ncpu)(delayed(model)(x, -1) for x in m)
    # r1 = Parallel(n_jobs=ncpu)(delayed(model)(x, 1) for x in m)
    Timer(t1)
    # np.savez(datafile, r0 = r0, r1 = r1)
    np.savez(datafile, r0 = r0, m = m)

# Load tau
r = np.load(datafile)
f_tau = r['r0']
m_tau = r['m']

# Load OmB
r0 = np.load('data/1.npz')
m0 = r0['m']
f0 = r0['f']

# Load T21
T21 = np.load('data/4.npz')
m_t21 = T21['m']
f_t21 = T21['f']
f_t21_ = T21['f1']

# Load Tk
Tk = np.load('data/7.npz')
m_tk = Tk['m']
f_tk = Tk['r0']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

plt.loglog(m0, f0, 'k', linewidth=LineWidth, label = '$\delta \Omega_{\mathrm{b}}/\Omega_{\mathrm{b}} < 30\%$')
plt.loglog(m_tau, f_tau, 'r', linewidth=LineWidth, label = '$\\tau_{\mathrm{rei}} < 0.0561$')
plt.loglog(m_t21, f_t21, 'g', linewidth=LineWidth, label = '$T_{21} < -100{\mathrm{mK}}$')
plt.loglog(m_t21, f_t21_, '--g', linewidth=LineWidth, label = '$T_{21} < -100{\mathrm{mK}}$, no Evo')
plt.loglog(m_tk, f_tk, 'b', linewidth=LineWidth, label = 'LyA $T_{\mathrm{k}}$ Limits')

plt.xlabel('$m_0\ [m_{\odot}]$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$f_{\mathrm{bh}}$',fontsize=FontSize,fontname='Times New Roman')
plt.legend(fontsize=FontSize/1.5,loc = 'lower left')
plt.title('Monochromatic PBH',fontsize=FontSize)

plt.xlim([1e-1, 1e4])
plt.ylim([1e-10, 1.5])

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png',bbox_inches='tight',dpi=1000)
