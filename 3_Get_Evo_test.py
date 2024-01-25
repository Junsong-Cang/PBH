type = 2

reload = 0
datafile = 'data/Get_Evo_test.npz'
LineWidth = 2
FontSize = 15

from src.main import *
from matplotlib.colors import LogNorm

if reload:
    t1 = TimeNow()
    r = Get_Evo(show_status = 1)
    Timer(t1)
    np.savez(datafile, MF = r['MF'], m_vec = r['m_vec'], z_vec = r['z_vec'], 
             omb_ratio = r['omb_ratio'], m0_vec = r['m0_vec'], fbh_ratio = r['fbh_ratio'], m_ratio = r['m_ratio'])
    
r = np.load(datafile)
mf = r['MF']
m_vec = r['m_vec']
z = r['z_vec']
b = r['omb_ratio']
m0 = r['m0_vec']
fbh_ratio = r['fbh_ratio']
m_ratio = r['m_ratio']
nz = len(z)

# Get plot
plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

if type == 1:
    
    # 2D MF

    # get the axis right
    PlotData = np.empty(np.shape(mf))
    for zid in np.arange(0,nz):
        m_axis = m_vec[zid,:]
        phi_axis = mf[zid,:]
        PlotData[zid,:] = np.interp(np.log(m0), np.log(m_axis), phi_axis, left = np.min(mf), right = np.max(mf))[:]
    
    m0, z = np.meshgrid(m0, z)
    fig,ax = plt.subplots()
    c=ax.pcolor(m0, z, PlotData, cmap='jet',norm = LogNorm(vmin=PlotData.min(), vmax=PlotData.max()))
    # c=ax.pcolor(m0, z, PlotData, cmap='jet')
    
    plt.xlabel('$m\ [m_{\odot}]$',fontsize=FontSize,fontname='Times New Roman')
    plt.ylabel('$z$',fontsize=FontSize,fontname='Times New Roman')
    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(size=FontSize)
    plt.yticks(size=FontSize)
    clb = plt.colorbar(c)
    plt.title('Mass Function $\phi \equiv {\mathrm{d}}{\mathrm{ln}}f_{\mathrm{bh}}/{\mathrm{d}}{\mathrm{ln}} m$',fontsize=FontSize)
    clb.ax.set_title('$\Phi$')
    # plt.title(TitleStr[PlotType],fontsize=FontSize)
    plt.tight_layout()
    plt.savefig('/Users/cangtao/Desktop/tmp.png',dpi=500)

elif type == 2:
    # 1D Phi shape
    Normalise = 1
    z1 = 600
    z2 = 25
    z3 = 10
    z4 = 8
    zid1 = 0
    zid2 = Find_Index(z2, z)
    zid3 = Find_Index(z3, z)
    zid4 = nz-1
    m1 = m_vec[zid1,:]
    p1 = mf[zid1,:]
    m2 = m_vec[zid2,:]
    p2 = mf[zid2,:]
    m3 = m_vec[zid3,:]
    p3 = mf[zid3,:]
    m4 = m_vec[zid4,:]
    p4 = mf[zid4,:]
    
    normalise = lambda x: x/np.max(x)

    if Normalise:
        p1 = normalise(p1)
        p2 = normalise(p2)
        p3 = normalise(p3)
        p4 = normalise(p4)
    # print(m2)
    plt.plot(m1,  p1, 'k', linewidth=LineWidth, label = 'Initial')
    plt.plot(m2,  p2, 'r', linewidth=LineWidth, label = '$z = 25$')
    plt.plot(m3,  p3, 'b', linewidth=LineWidth, label = '$z = 10$')
    plt.plot(m4,  p4, 'g', linewidth=LineWidth, label = '$z = 8$')

    plt.xlabel('$m\ [m_{\odot}]$',fontsize=FontSize,fontname='Times New Roman')
    plt.ylabel('$\Phi$',fontsize=FontSize,fontname='Times New Roman')
    plt.legend(fontsize=FontSize,loc = 'lower right')
    plt.xscale('log')
    if not Normalise:
        plt.yscale('log')
    
    # plt.title('re-scaled',fontsize=FontSize)
    plt.xticks(size=FontSize)
    plt.yticks(size=FontSize)
    plt.tight_layout()

    plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=500)

elif type == 3:
    # OmB
    plt.plot(z,  b, 'k', linewidth=LineWidth, label = '$\Omega_{\mathrm{b}}$')
    plt.plot(z,  fbh_ratio, 'r', linewidth=LineWidth, label = '$f_{\mathrm{bh}}$')
    plt.xlabel('$z$',fontsize=FontSize,fontname='Times New Roman')
    plt.ylabel('ratio',fontsize=FontSize,fontname='Times New Roman')
    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(size=FontSize)
    plt.yticks(size=FontSize)
    plt.legend(fontsize=FontSize,loc = 'upper right')
    plt.tight_layout()
    plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=500)
    print(np.min(b))

    

