reload = 0
EoR = 1
m1 = 1e0
m2 = 1e4
nm = 100

from src.Radio import *
m = np.logspace(np.log10(m1), np.log10(m2), nm)

if reload:
    t1 = TimeNow()
    f = np.zeros(nm)
    for idx in np.arange(0, nm):
        f[idx] = Find_fbh_Radio(mc = m[idx], sbh = -1, Use_EoR = EoR)
        print(idx/nm)
    Timer(t1)
    np.savez('tmp.npz', f = f)

f = np.load('tmp.npz')['f']

plt.loglog(m,f)
plt.show()
