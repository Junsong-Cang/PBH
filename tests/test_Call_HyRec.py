from src.Recombination import *

reload = 1
m = 1e-3
LineWidth = 1.4
FontSize = 12

if reload:
    r0 = Call_HyRec(fbh0 = 0, mc = m, sbh = -1, show_status = 0)
    r1 = Call_HyRec(fbh0 = 0.99, mc = m, sbh = -1, show_status = 1)

    z0 = r0['z']
    x0 = r0['xe']
    t0 = r0['t']
    b0 = r0['OmB_ratio']
    zb0 = r0['zb']
    
    z1 = r1['z']
    x1 = r1['xe']
    t1 = r1['t']
    b1 = r1['OmB_ratio']
    zb1 = r1['zb']
    
    np.savez('tmp.npz', z0 = z0, x0 = x0, t0 = t0, b0 = b0, zb0 = zb0, z1 = z1, x1 = x1, t1 = t1, b1 = b1, zb1 = zb1)

r = np.load('tmp.npz')
z0 = r['z0']
x0 = r['x0']
t0 = r['t0']
b0 = r['b0']
zb0 = r['zb0']

z1 = r['z1']
x1 = r['x1']
t1 = r['t1']
b1 = r['b1']
zb1 = r['zb1']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, axs = plt.subplots(1, 3, sharex = False, sharey = False)
fig.set_size_inches(10, 3)

axs[0].loglog(1+z0, t0, 'k', linewidth = LineWidth, label='LCDM')
axs[0].loglog(1+z1, t1, 'r', linewidth = LineWidth, label='PBH')
axs[0].set_title('$T_{\mathrm{k}}$ (K)',fontsize=FontSize)
axs[0].legend(fontsize=FontSize, loc = 'lower right')
axs[0].set_xlabel('$1+z$',fontsize=FontSize,fontname='Times New Roman')
axs[0].tick_params(axis='both', which='both', labelsize = FontSize)
#axs[0].set_xlim(-3.2, 3.2)
#axs[0].set_ylim(-1.2, 1.2)

axs[1].loglog(1+z0, x0, 'k', linewidth = LineWidth)
axs[1].loglog(1+z1, x1, 'r', linewidth = LineWidth)
axs[1].set_title('$x_{\mathrm{e}}$',fontsize=FontSize)
axs[1].set_xlabel('$1+z$',fontsize=FontSize,fontname='Times New Roman')
axs[1].tick_params(axis='both', which='both', labelsize = FontSize)
#axs[1].set_xlim(-3.2, 3.2)
#axs[1].set_ylim(0, 1.2)

axs[2].loglog(1+zb0, b0, 'k', linewidth = LineWidth)
axs[2].loglog(1+zb1, b1, 'r', linewidth = LineWidth)
axs[2].set_title('$\Omega_{\mathrm{b}} / \Omega_{\mathrm{b0}}$',fontsize=FontSize)
axs[2].set_xlabel('$1+z$',fontsize=FontSize,fontname='Times New Roman')
axs[2].tick_params(axis='both', which='both', labelsize = FontSize)
#axs[1].set_xlim(-3.2, 3.2)
#axs[1].set_ylim(0, 1.2)

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)
t17 = np.interp(x = 17, xp = z1[::-1], fp = t1[::-1])
print(t17)
