from src.mdot_kernel import *

n = 20

from PyLab import *

t1 = TimeNow()
for idx in np.arange(0, n):
    r = mdot_kernel()
    print(idx)
    Timer(t1)

