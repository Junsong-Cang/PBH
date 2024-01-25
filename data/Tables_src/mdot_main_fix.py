nm = 300
nz = 150
z1 = 6.99
z2 = 1000.1
m1 = 1e-3
m2 = 1e5

reload = 1
ncpu = 12
EoR_z = 36

swap_file_old = 'swap/mdot_tab_main_swap.npz'
swap_file = 'swap/mdot_tab_main_swap_fix.npz'

log_file = 'tmp.txt'

# Average speed:
# 1s per call

from src.mdot_kernel import *
N = nm*nz
indexs = np.arange(0, N)

lm1 = np.log10(m1)
lm2 = np.log10(m2)
lz1 = np.log10(z1)
lz2 = np.log10(z2)

m_vec = np.logspace(lm1, lm2, nm)
z_vec = np.logspace(lz1, lz2, nz)

params = np.empty((2, N))

idx = 0
for zid in np.arange(0, nz):
    for mid in np.arange(0, nm):
        params[0, idx] = z_vec[zid]
        params[1, idx] = m_vec[mid]
        idx = idx + 1

# loading prev results
r_old = np.load(swap_file_old)['r']
'''
mdot_vec = np.zeros((2, nz, nm))

idx = 0
for zid in np.arange(0,nz):
    for mid in np.arange(0,nm):
        mdot_vec[0, zid, mid] = r_old[idx, 0]
        mdot_vec[1, zid, mid] = r_old[idx, 1]
        idx_ = round(r_old[idx, 2])
        if not (idx == idx_):
            raise Exception
        idx = idx + 1
'''

def model(id = 1):

    z = params[0, id]
    m = params[1, id]
    r = []
    # r1 = mdot_kernel(z = z, M_PBH =m, Use_EoR = False)
    r1 = r_old[id,0]

    if z > EoR_z:
        r2 = r1
    else:
        r2 = mdot_kernel(z = z, M_PBH =m, Use_EoR = True)
        SaySomething(log_file)

    r.append(r1)
    r.append(r2)
    r.append(id)
    
    r = np.array(r)
    
    return r

if reload:
    t1 = TimeNow()
    r = Parallel(n_jobs=ncpu)(delayed(model)(x) for x in indexs)
    Timer(t1)
    np.savez(swap_file, r = r)

