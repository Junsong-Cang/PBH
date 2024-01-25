from src.main import *

LineWidth = 2
FontSize = 15

reload = 1

nz = 96
m = 0.023
z1 = 7
z2 = 1000
Use_EoR = 1
ncpu = 12
datafile = 'data/mdot_interp_1D.npz'

z = np.logspace(np.log10(z1), np.log10(z2), nz)

f0 = lambda x : Get_mdot(z = x, m = m, Use_Edd_Limit = 0, Use_EoR = Use_EoR, Do_not_Interpolate = True)
f1 = lambda x : Get_mdot(z = x, m = m, Use_Edd_Limit = 0, Use_EoR = Use_EoR)

if reload:
    t1 = TimeNow()    
    r0 = Parallel(n_jobs=ncpu)(delayed(f0)(x) for x in z)
    r1 = Parallel(n_jobs=ncpu)(delayed(f1)(x) for x in z)
    Timer(t1)
    np.savez(datafile, r0 = r0, r1 = r1)

r = np.load(datafile)
r0 = r['r0']
r1 = r['r1']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()

plt.plot(z, r0, 'k', linewidth=LineWidth, label = 'No Interpolation')
plt.plot(z, r1, '--r', linewidth=LineWidth, label = 'Interpolation')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$z$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$\dot{m}$',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.legend(fontsize=FontSize,loc = 'lower right')
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)
