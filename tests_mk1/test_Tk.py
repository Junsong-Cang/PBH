from PyLab import *

nz = 100
r1 = np.zeros(nz)
r2 = np.zeros(nz)

z = np.logspace(0, 3.2, nz) - 1

from src.mdot_kernel_mk1 import *

for idx in np.arange(0, nz):
    r1[idx] = Tk_ez(z[idx])
    xe, r2[idx] = LCDM_HyRec(z = z[idx])

plt.loglog(1+z, r1, 'k')
plt.loglog(1+z, r2, 'r')
plt.show()
