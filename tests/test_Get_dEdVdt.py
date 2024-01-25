from src.Recombination import *

LineWidth = 1.7
FontSize = 17

reload = 1

if reload:
    r = Get_dEdVdt()
    np.savez('tmp.npz', z = r['z'], d1 = r['dEdVdt_HIon'], d3 = r['dEdVdt_LyA'], d0 = r['dEdVdt_inj'],
             d4 = r['dEdVdt_Heat'], f1 = r['fc_HIon'], f3 = r['fc_LyA'], f4 = r['fc_Heat'])

r = np.load('tmp.npz')

z = r['z']
d0 = r['d0']
d1 = r['d1']
d3 = r['d3']
d4 = r['d4']

f1 = r['f1']
f3 = r['f3']
f4 = r['f4']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True

fig, axs = plt.subplots(1, 2, sharex = False, sharey = False)
fig.set_size_inches(10, 4)

axs[0].loglog(z, d0, 'k', linewidth = LineWidth, label='Injection')
axs[0].loglog(z, d1, 'r', linewidth = LineWidth, label='HIon')
axs[0].loglog(z, d3, 'g', linewidth = LineWidth, label='LyA')
axs[0].loglog(z, d4, 'b', linewidth = LineWidth, label='Heat')

axs[0].legend(fontsize=FontSize, loc = 'upper left')
axs[0].set_title('Deposition Rate',fontsize=FontSize)
axs[0].set_xlabel('$z$',fontsize=FontSize,fontname='Times New Roman')
axs[0].set_ylabel('d$E$/d$V$d$t\ [{\mathrm{eVcm^{-3}s^{-1}}}]$',fontsize=FontSize,fontname='Times New Roman')
axs[0].tick_params(axis='both', which='both', labelsize = FontSize)
#axs[0].set_xlim(-3.2, 3.2)
#axs[0].set_ylim(-1.2, 1.2)

axs[1].loglog(z, f1, 'r', linewidth = LineWidth)
axs[1].loglog(z, f3, 'g', linewidth = LineWidth)
axs[1].loglog(z, f4, 'b', linewidth = LineWidth)

axs[1].set_title('Depositon Efficiency',fontsize=FontSize)
axs[1].set_xlabel('$z$',fontsize=FontSize,fontname='Times New Roman')
axs[1].set_ylabel('$f_{\mathrm{c}}$',fontsize=FontSize,fontname='Times New Roman')

#axs[1].set_xlim(-3.2, 3.2)
#axs[1].set_ylim(0, 1.2)

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)

print(f1)