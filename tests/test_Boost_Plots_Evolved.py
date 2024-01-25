from src.main import *

reload = 0
m1 = 1e0
m2 = 1e1
m3 = 1e3
nz = 2000

LineWidth = 2
FontSize = 25

EoR = 1
nz0 = 100
z0 = np.linspace(7, 40, nz0)

if reload:
    b = np.zeros((3, nz))
    m0 = [m1, m2, m3]
    Evo_1 = Get_Evo_Lite(fbh0 = 1e-10, Use_EoR = EoR, nz = nz, m0 = m1)
    Evo_2 = Get_Evo_Lite(fbh0 = 1e-10, Use_EoR = EoR, nz = nz, m0 = m2)
    Evo_3 = Get_Evo_Lite(fbh0 = 1e-10, Use_EoR = EoR, nz = nz, m0 = m3)
    Evos = [Evo_1, Evo_2, Evo_3]
    
    for idx in np.arange(0, 3):
        Evo = Evos[idx]
        z = Evo['z']
        m = m0[idx] * Evo['m_ratio']
        print(idx)
        
        for zid in np.arange(nz):
            z_ = z[zid]
            m_ = m[zid]
            if z_ < 40:
                b[idx, zid] = BoostFactor(z = z_, m = m_, Use_EoR = EoR)
            else:
                b[idx, zid] = 1
    
    b0 = np.zeros(nz0)
    for zid in np.arange(nz0):
        print(zid)
        z_ = z0[zid]
        b0[zid] = Boost_Factor_Lite(z = z_, Use_EoR = EoR, ncpu = 11)
        
    np.savez('tmp.npz', z = z, b = b, b0 = b0)

r = np.load('tmp.npz')

b0 = r['b0']
z = r['z']

b = r['b']
b1 = b[0,:]
b2 = b[1,:]
b3 = b[2,:]

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

plt.plot(z0, b0, 'k', linewidth=LineWidth, label = 'No Edd')
plt.plot(z, b1, 'r', linewidth=LineWidth, label = '$m_0 = 1$')
plt.plot(z, b2, 'g', linewidth=LineWidth, label = '$m_0 = 10^{1}$')
plt.plot(z, b3, 'b', linewidth=LineWidth, label = '$m_0 = 10^{3}$')
plt.yscale('log')
plt.xlabel('$z$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$B$',fontsize=FontSize,fontname='Times New Roman')
plt.legend(fontsize=FontSize,loc = 'upper right')
#plt.title('$\delta \Omega_{\mathrm{b}}<30 \%$',fontsize=FontSize)

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.xlim([7, 40])

plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png',bbox_inches='tight',dpi=1000)

print(b3)
