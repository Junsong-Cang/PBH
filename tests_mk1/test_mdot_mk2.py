veff_model = 5
m = 100000
nz = 50

from src.mdot_kernel_mk2 import *

DataPath = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/'
if m == 100000:
    f = DataPath + '2003.02778.fig1.m1e5.txt'
elif m == 10000:
    f = DataPath + '2003.02778.fig1.m1e4.txt'
elif m == 300:
    f = DataPath + '2003.02778.fig1.m3e2.txt'
elif m == 100:
    f = DataPath + '2003.02778.fig1.m1e2.txt'
elif m == 10:
    f = DataPath + '2003.02778.fig1.m10.txt'
elif m == 1:
    f = DataPath + '2003.02778.fig1.m1.txt'

zp, r1 = Read_Curve(
        File = f,
        nx = nz,
        model = 2,
        Convert_x = 1,
        Convert_y = 1)

r2 = np.zeros(nz)
for veff_model in np.arange(0, 5):
    for idx in np.arange(0, nz):
        z = zp[idx]-1
        r2[idx] = mdot_kernel(m = m, z = z, veff_model = veff_model)
    plt.loglog(zp, r2, label = str(veff_model))
    

plt.loglog(zp, r1, 'k', label = 'Paper')
#plt.legend(fontsize=FontSize,loc = 'lower left')
plt.legend(loc = 'lower right')

print(r2)
plt.show()
