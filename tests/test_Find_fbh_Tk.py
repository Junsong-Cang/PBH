from src.Heating_Limits import *

reload = 1
m1 = 1
m2 = 300
nm = 100

m = np.logspace(np.log10(m1), np.log10(m2), nm)

if reload:
    t1 = TimeNow()
    r = Parallel(n_jobs=12)(delayed(Find_fbh_Tk)(mc = x, sbh = -1) for x in m)
    Timer(t1)
    np.savez('tmp.npz', r = r)

r = np.load('tmp.npz')['r']

plt.loglog(m,r)
plt.show()
