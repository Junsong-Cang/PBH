'''
from PyLab import *

z = np.logspace(0, 3.5, 1000) - 1
x0, t0 = LCDM_HyRec(z=z, Use_EoR=0)
x1, t1 = LCDM_HyRec(z = z, Use_EoR=1)

plt.loglog(z, x0, 'k')
plt.loglog(z, x1, '--r')
plt.show()
'''

'''
x = 1
while True:
    x = x+1
    print('--', x)
    if x > 100:
        break
    print(x)

'''

n = 200
nr = 100
nz = 100

N = n**3 * nr * nz
print(N)
