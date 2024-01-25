from src.main import *

reload = 1
m1 = 1e-3
m2 = 8

z = np.logspace(np.log10(8), 3, 100)
LineWidth = 2
FontSize = 18

if reload:
    t1 = TimeNow()
    r1 = Parallel(n_jobs=12)(delayed(Get_mdot)(z = x, m = m1, Use_EoR = 0) for x in z)
    r2 = Parallel(n_jobs=12)(delayed(Get_mdot)(z = x, m = m1, Use_EoR = 1) for x in z)
    r3 = Parallel(n_jobs=12)(delayed(Get_mdot)(z = x, m = m2, Use_EoR = 0) for x in z)
    r4 = Parallel(n_jobs=12)(delayed(Get_mdot)(z = x, m = m2, Use_EoR = 2) for x in z)
    '''
    r1_ = Parallel(n_jobs=12)(delayed(Get_mdot)(z = x, m = m1, Use_EoR = 0, Do_not_Interpolate = 1) for x in z)
    r2_ = Parallel(n_jobs=12)(delayed(Get_mdot)(z = x, m = m1, Use_EoR = 1, Do_not_Interpolate = 1) for x in z)
    r3_ = Parallel(n_jobs=12)(delayed(Get_mdot)(z = x, m = m2, Use_EoR = 0, Do_not_Interpolate = 1) for x in z)
    r4_ = Parallel(n_jobs=12)(delayed(Get_mdot)(z = x, m = m2, Use_EoR = 1, Do_not_Interpolate = 1) for x in z)
    '''

    Timer(t1)
    # np.savez('tmp.npz', r1 = r1, r2 = r2, r3 = r3, r4 = r4, r1_ = r1_, r2_ = r2_, r3_ = r3_, r4_ = r4_)
    np.savez('tmp.npz', r1 = r1, r2 = r2, r3 = r3, r4 = r4)


d = np.load('tmp.npz')
r1, r2, r3, r4 = d['r1'], d['r2'], d['r3'], d['r4']
# r1_, r2_, r3_, r4_ = d['r1_'], d['r2_'], d['r3_'], d['r4_']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots()


plt.plot(z, r2, 'k', linewidth=LineWidth, label = 'm1, EoR')
plt.plot(z, r1, '--k', linewidth=LineWidth, label = 'm1, no EoR')
plt.plot(z, r4, 'b', linewidth=LineWidth, label = 'm2, EoR')
plt.plot(z, r3, '--b', linewidth=LineWidth, label = 'm2, no EoR')

'''
plt.plot(z, r1_, '+k', linewidth=LineWidth)
plt.plot(z, r2_, '+r', linewidth=LineWidth)
plt.plot(z, r3_, '+b', linewidth=LineWidth)
plt.plot(z, r4_, '+g', linewidth=LineWidth)
'''

plt.xscale('log')
plt.yscale('log')

plt.xlabel('$z$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$y$',fontsize=FontSize,fontname='Times New Roman')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.legend(fontsize=FontSize,loc = 'upper right')
plt.tight_layout()

plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=500)
