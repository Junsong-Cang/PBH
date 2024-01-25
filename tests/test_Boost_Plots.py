from src.main import *

reload = 1
m1 = 1e0
m2 = 1e1
m3 = 1e3

LineWidth = 2
FontSize = 25

nz = 100
z = np.linspace(8, 39, nz)
EoR = 1

if reload:
    b0 = np.zeros(nz)
    b1 = np.zeros(nz)
    b2 = np.zeros(nz)
    b3 = np.zeros(nz)
    t1 = TimeNow()
    r = np.load('tmp.npz')
    b0 = r['b0']

    for idx in np.arange(0, nz):
        z_ = z[idx]
        #b0[idx] = Boost_Factor_Lite(z = z_, Use_EoR = EoR)
        b1[idx] = BoostFactor(z = z_, m = m1, Use_EoR = EoR)
        b2[idx] = BoostFactor(z = z_, m = m2, Use_EoR = EoR)
        b3[idx] = BoostFactor(z = z_, m = m3, Use_EoR = EoR)
        print(idx/nz)
    Timer(t1)
    np.savez('tmp.npz', b0 = b0, b1 = b1, b2 = b2, b3 = b3)

r = np.load('tmp.npz')

b0 = r['b0']
b1 = r['b1']
b2 = r['b2']
b3 = r['b3']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

plt.plot(z, b0, 'k', linewidth=LineWidth, label = 'No Edd')
plt.plot(z, b1, 'r', linewidth=LineWidth, label = '$m = 1$')
plt.plot(z, b2, 'g', linewidth=LineWidth, label = '$m = 10^{1}$')
plt.plot(z, b3, 'b', linewidth=LineWidth, label = '$m = 10^{3}$')
plt.yscale('log')
plt.xlabel('$z$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$B$',fontsize=FontSize,fontname='Times New Roman')
plt.legend(fontsize=FontSize,loc = 'upper right')
#plt.title('$\delta \Omega_{\mathrm{b}}<30 \%$',fontsize=FontSize)

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png',bbox_inches='tight',dpi=1000)

print(b3)
