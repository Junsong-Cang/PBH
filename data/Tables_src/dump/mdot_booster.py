nm = 100
nz = 100
m1 = 1e4
m2 = 1e5
z1 = 7
z2 = 1000
ncpu = 12
reload = 1
swap_file = 'mdot_booster.npz'
log_file = 'mdot_booster.txt'

from src.src_1 import *
N = nm*nz
lm1 = np.log10(m1)
lm2 = np.log10(m2)
lz1 = np.log10(z1)
lz2 = np.log10(z2)

m_vec = np.logspace(lm1, lm2, nm)
z_vec = np.logspace(lz1, lz2, nz)

t1 = TimeNow()

model = lambda x : Get_mdot(z = x, m = 1e10, Use_Edd_Limit = 0, Do_not_Interpolate=1)
r = Parallel(n_jobs=ncpu)(delayed(model)(x) for x in z_vec)
Timer(t1)    

print(np.min(r))
plt.loglog(z_vec,r)
plt.show()
