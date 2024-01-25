from src.mdot_kernel import *

reload = 1
nz = 50

from PyLab import *
z = np.logspace(0, 2.1, nz)

def model(z_):
    r = mdot_kernel(
        z = z_,
        M_PBH = 1E1,
        Use_EoR = 0
    )
    # print('----')
    return r

if reload:
    t1 = TimeNow()
    # r = Parallel(n_jobs=11)(delayed(model)(z_) for z_ in z)
    r = np.zeros(nz)
    for idx in np.arange(0, nz):
        print(idx)
        r[idx] = model(z_ = z[idx])
        Timer(t1)
    np.savez('tmp.npz', r = r)

'''
r = np.load('tmp.npz')['r']
for idx in np.arange(0, nz):
    if r[idx]>10:
        r[idx] = 10
plt.loglog(1+z, r)
plt.show()
'''
