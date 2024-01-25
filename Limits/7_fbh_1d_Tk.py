from src.Heating_Limits import *

reload = 0
m1 = 1
m2 = 300
nm = 50
datafile = 'data/7.npz'

m = np.logspace(np.log10(m1), np.log10(m2), nm)
def model(x, s):
    if s<0:
        r = Find_fbh_Tk(mc = x, sbh = s, zmax = 1000)
    else:
        r = Find_fbh_Tk(mc = x, sbh = s)
    SaySomething('tmp.txt')
    return r

if reload:
    t1 = TimeNow()
    r0 = Parallel(n_jobs=11)(delayed(model)(x, -1) for x in m)
    # r1 = Parallel(n_jobs=11)(delayed(model)(x, 1) for x in m)
    Timer(t1)
    # np.savez(datafile, r0 = r0, r1 = r1)
    np.savez(datafile, r0 = r0, m = m)

r = np.load(datafile)
r0 = r['r0']
#r1 = r['r1']

plt.loglog(m, r0, 'k')
#plt.semilogx(m, r1, 'r')
plt.show()
