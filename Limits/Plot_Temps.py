from src.Recombination import *
from src.Radio import *

fbh0 = 1e-2
#z = np.logspace(0, 3, 1000) - 1
m1 = 1e1
m2 = 1e3

LineWidth = 2
FontSize = 18

# --------
Rec1 = Call_HyRec(fbh0 = fbh0, sbh = -1, mc = m1, zmax = 1000, zmin = 8, Use_EoR = False)
Rec2 = Call_HyRec(fbh0 = fbh0, sbh = -1, mc = m2, zmax = 1000, zmin = 8, Use_EoR = False)
R1 = T_Radio(mc = m1, sbh = -1, fbh0 = fbh0, Use_EoR = True)
R2 = T_Radio(mc = m2, sbh = -1, fbh0 = fbh0, Use_EoR = True)

z1, x1, t1 = Rec1['z'], Rec1['xe'], Rec1['t']
z2, x2, t2 = Rec2['z'], Rec2['xe'], Rec2['t']
zr1, tr1 = R1['z'], R1['T']
zr2, tr2 = R2['z'], R2['T']

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, axs = plt.subplots(1, 2, sharex = False, sharey = False)
fig.set_size_inches(10, 4)

#axs[0].loglog(1+z1, x1, 'k', linewidth = LineWidth, label='$T{\mathrm{k}}$')
axs[0].loglog(1+z1, x1, 'k', linewidth = LineWidth, label='$m_0 = 10$')
axs[0].loglog(1+z2, x2, 'r', linewidth = LineWidth, label='$m_0 = 10^3$')

axs[0].legend(fontsize=FontSize, loc = 'lower left')
axs[0].set_title('$x_{\mathrm{e}}$',fontsize=FontSize)
axs[0].set_xlabel('$1 + z$',fontsize=FontSize,fontname='Times New Roman')
#axs[0].set_ylabel('$y$',fontsize=FontSize,fontname='Times New Roman')
axs[0].tick_params(axis='both', which='both', labelsize = FontSize)
# Set axis limits
axs[0].set_xlim(1, 2000)
#axs[0].set_ylim(-1.2, 1.2)

axs[1].loglog(1+z1, t1, 'k', linewidth = LineWidth, label='$T_{\mathrm{k}}$')
axs[1].loglog(1+z2, t2, 'r', linewidth = LineWidth)
axs[1].loglog(1+zr1, tr1, '--k', linewidth = LineWidth, label='$T_{\mathrm{Radio}}$')
axs[1].loglog(1+zr2, tr2, '--r', linewidth = LineWidth)
axs[1].errorbar([5.3, 5.8, 7.08], [10000, 7980, 10**4.21], [1000, 2500, 2381],color = 'k', fmt='+', linewidth=LineWidth)

axs[1].legend(fontsize=FontSize, loc = 'upper right')

axs[1].set_title('Temperatures',fontsize=FontSize)
axs[1].set_xlabel('$1+z$',fontsize=FontSize,fontname='Times New Roman')
axs[1].set_xlim(1, 1000)

#axs[1].set_xlim(-3.2, 3.2)
#axs[1].set_ylim(0, 1.2)

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.png', dpi=1000)

print(z1[0])
z1 = z1[::-1]
z2 = z2[::-1]
x1 = x1[::-1]
x2 = x2[::-1]


t1 = Find_Optical_Depth(z = z1, xe = x1, zmax = 100)
t2 = Find_Optical_Depth(z = z2, xe = x2, zmax = 100)

print(t1, t2)
plt.show()