from src.src_1 import *

reload = 1
PlotType = 2
f0 = 1e-5
nm = 500
nz = 1000
datafile = 'data/test_1.npz'
zmax = 1000
zmin = 8
ncpu = 12
m1 = 1
m2 = 100

LineWidth = 2
FontSize = 15

from matplotlib.colors import LogNorm

m0 = np.logspace(np.log10(m1), np.log10(m2), nm)
# m0 = np.linspace(m1, m2, nm)

fun = lambda x, usb, Use_EoR: Get_Evo_Lite(fbh0 = f0, m0 = x, zmin = zmin, zmax = zmax, nz = nz, Use_OmB_FeedBack = usb, Use_EoR = Use_EoR)

if reload:
    t1 = TimeNow()

    d1 = Parallel(n_jobs = ncpu)(delayed(fun)(x, False, False) for x in m0)
    d2 = Parallel(n_jobs = ncpu)(delayed(fun)(x, False, True) for x in m0)
    d3 = Parallel(n_jobs = ncpu)(delayed(fun)(x, True, False) for x in m0)
    d4 = Parallel(n_jobs = ncpu)(delayed(fun)(x, True, True) for x in m0)

    Timer(t1)
    
    # Organize
    z = d1[0]['z']
    '''
    for mid in np.arange(0, nm):
        m1[:,mid] = d1[mid]['m'][:]
        m2[:,mid] = d2[mid]['m'][:]
        b[:,mid] = d1[mid]['OmB'][:]
    '''
    # np.savez(datafile, z = z, m1 = m1, m2 = m2, b = b)
    np.savez(datafile, z = z, d1 = d1, d2 = d2, d3 = d3, d4 = d4)

d = np.load(datafile, allow_pickle = True)
d1 ,d2, d3, d4= d['d1'], d['d2'], d['d3'], d['d4']

m1 = np.empty((nz, nm))
m2 = np.empty((nz, nm))
m3 = np.empty((nz, nm))
m4 = np.empty((nz, nm))
b3 = np.empty((nz, nm))
b4 = np.empty((nz, nm))

z = d1[0]['z']
for mid in np.arange(0, nm):
    m1[:,mid] = d1[mid]['m_ratio'][:]
    m2[:,mid] = d2[mid]['m_ratio'][:]
    m3[:,mid] = d3[mid]['m_ratio'][:]
    m4[:,mid] = d4[mid]['m_ratio'][:]
    b3[:,mid] = d3[mid]['OmB_ratio'][:]
    b4[:,mid] = d3[mid]['OmB_ratio'][:]

PlotDataSet = [m1, m2, m3, m4, b3, b4]
PlotData = PlotDataSet[PlotType]
cb_str = ['$m/m_0$', '$m/m_0$', '$m/m_0$', '$m/m_0$', '$\Omega_{\mathrm{b}}/\Omega_{\mathrm{b0}}$', '$\Omega_{\mathrm{b}}/\Omega_{\mathrm{b0}}$']
TitleStr = ['No OmB Feedback, no EoR', 'No OmB Feedback', 'No EoR', 'Full', 'No EoR', 'Full']

# Get plot
plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

m0, z = np.meshgrid(m0, z)
fig,ax = plt.subplots()
if PlotType<4:
    c=ax.pcolor(m0, z, PlotData, cmap='jet',norm = LogNorm(vmin=PlotData.min(), vmax=PlotData.max()))
else:
    c=ax.pcolor(m0, z, PlotData, cmap='jet')
    
plt.xlabel('$m_0\ [m_{\odot}]$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$z$',fontsize=FontSize,fontname='Times New Roman')
plt.xscale('log')
plt.yscale('log')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
clb = plt.colorbar(c)
clb.ax.set_title(cb_str[PlotType])
plt.title(TitleStr[PlotType],fontsize=FontSize)
plt.tight_layout()
# plt.savefig('/Users/cangtao/Desktop/tmp.png',dpi=500)
plt.show()
