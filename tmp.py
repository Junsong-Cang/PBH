reload = 1
LineWidth = 2
FontSize = 18

from p21c_tools import *

z = np.linspace(8, 20, 200)
if reload:
    t1 = TimeNow()
    f0 = fcoll_function(z = z, Use_EoR = 0)
    f1 = fcoll_function(z = z, Use_EoR = 1)
    Timer(t1)
    np.savez('tmp.npz', f0 = f0, f1 = f1)

r = np.load('tmp.npz')
f0 = r['f0']
f1 = r['f1']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

plt.plot(z, f0, 'k', linewidth=LineWidth, label = 'no EoR')
plt.plot(z, f1, 'r', linewidth=LineWidth, label = 'with EoR')
plt.yscale('log')

plt.xlabel('$z$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$f_{\mathrm{coll}}$',fontsize=FontSize,fontname='Times New Roman')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.title('Fraction of matter collapsed into halos',fontsize=FontSize)
plt.legend(fontsize=FontSize,loc = 'lower left')
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)
