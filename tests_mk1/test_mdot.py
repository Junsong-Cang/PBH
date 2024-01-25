mid = eval(input('enter mid:'))
nz = 100
xe = 1

from src.mdot_kernel_mk1 import *
m_vec = [100000, 10000, 1000, 300]
m = m_vec[mid]
z1 = 2
z2 = 5000
z = np.logspace(np.log10(1+z1), np.log10(1+z2), nz)
DataPath = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/'
if m == 100000 and xe < 1e-2:
    f = DataPath + '0709.0524.Fig4a.top1.txt'
elif m == 10000 and xe < 1e-2:
    f = DataPath + '0709.0524.Fig4a.top2.txt'
elif m == 1000 and xe < 1e-2:
    f = DataPath + '0709.0524.Fig4a.top3.txt'
elif m == 300 and xe < 1e-2:
    f = DataPath + '0709.0524.Fig4a.top4.txt'
elif m == 100000 and xe > 1e-2:
    f = DataPath + '0709.0524.Fig4b.top1.txt'
elif m == 10000 and xe > 1e-2:
    f = DataPath + '0709.0524.Fig4b.top2.txt'
elif m == 1000 and xe > 1e-2:
    f = DataPath + '0709.0524.Fig4b.top3.txt'
elif m == 300 and xe > 1e-2:
    f = DataPath + '0709.0524.Fig4b.top4.txt'

r1 = np.zeros(nz)
c = ['r', 'b', 'g', 'm', 'c', 'y']
for veff_model in np.arange(0, 6):
    for idx in np.arange(0, nz):
        r1[idx] = mdot_kernel(m = m, z = z[idx], xe = xe, veff_model = veff_model)
    plt.loglog(1+z, r1, color = c[veff_model], label = str(veff_model))
    
zp, r2 = Read_Curve(
        File = f,
        nx = 100,
        model = 2,
        Convert_x = 1,
        Convert_y = 1
        )
        
# plt.loglog(1+z, r1, 'k')
plt.loglog(zp, r2, 'k', label = 'Paper')
plt.legend(loc = 'lower right')

plt.show()
