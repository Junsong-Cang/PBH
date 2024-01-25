nm = 5
nz = 5
z1 = 1000
z2 = 2000
m1 = 1e-3
m2 = 1e5

reload = 1
ncpu = 12
EoR_z = 36
swap_file = 'mdot_tab_HiZ_swap.npz'
log_file = 'mdot_HiZ.txt'

from src.src_1 import *
N = nm*nz
lm1 = np.log10(m1)
lm2 = np.log10(m2)
lz1 = np.log10(z1)
lz2 = np.log10(z2)

m_vec = np.logspace(lm1, lm2, nm)
z_vec = np.logspace(lz1, lz2, nz)

indexs = np.arange(0, N)
params = np.empty((2, N))

idx = 0
for zid in np.arange(0, nz):
    for mid in np.arange(0, nm):
        params[0, idx] = z_vec[zid]
        params[1, idx] = m_vec[mid]
        idx = idx + 1

def model(id = 1):

    z = params[0, id]
    m = params[1, id]
    
    r = []
    r1 = Get_mdot(z = z, m = m, Do_not_Interpolate = True, Use_Edd_Limit = False, Use_EoR = False)
    r.append(r1)
    r.append(id)
    
    r = np.array(r)
    
    SaySomething(log_file)
    
    return r

if reload:
    t1 = TimeNow()
    r = Parallel(n_jobs=ncpu)(delayed(model)(x) for x in indexs)
    Timer(t1)
    np.savez(swap_file, r = r)
