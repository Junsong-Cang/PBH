from src.main import *

reload = 1
nm = 200
nz = 40
nb = 20

ncpu = 12

m1 = 1e-3
m2 = 1e4
z1 = 8
z2 = 40
bf1 = 1e-2
bf2 = 1.01
swap_file = '/Users/cangtao/cloud/GitHub/PBH/data/Boost_Tab_swap.npz'
data_file = '/Users/cangtao/cloud/GitHub/PBH/data/Boost_Tab.npz'

# Initialising
N = nm*nz*nb
print('Estimated time: ', N/3600)

OmB_LCDM = cosmo['OmB']
m_vec = np.logspace(np.log10(m1), np.log10(m2), nm)
bf_vec = np.logspace(np.log10(bf1), np.log10(bf2), nb)
z_vec = np.linspace(z1, z2, nz)
OmB_vec = OmB_LCDM*bf_vec

idx_vec = np.arange(0, N)

# I guess using a string will be more efficient
id = 0
params = np.empty((3, N))

for zid in np.arange(0, nz):
    for mid in np.arange(0, nm):
        for bid in np.arange(0, nb):
            z = z_vec[zid]
            m = m_vec[mid]
            OmB = OmB_vec[bid]
            params[0, id] = z
            params[1, id] = m
            params[2, id] = OmB
            id = id + 1

def model(idx):
    p = params[:,idx]
    z = p[0]
    m = p[1]
    OmB = p[2]
    r0 = BoostFactor(z = z, m = m, OmB = OmB, Use_EoR = 0, Use_Interp = False)
    r1 = BoostFactor(z = z, m = m, OmB = OmB, Use_EoR = 1, Use_Interp = False)
    r = np.array([r0, r1])
    SaySomething('tmp.txt')
    return r

if reload:
    t1 = TimeNow()
    r = Parallel(n_jobs=ncpu)(delayed(model)(x) for x in idx_vec)
    Timer(t1)
    np.savez(swap_file, r = r)

# Organising data
r = np.load(swap_file)['r']
print(np.shape(r))

# Main table, dimension : [EoR, OmB, z, m]

Boost_Tab = np.empty((2, nb, nz, nm))

idx = 0
for zid in np.arange(0, nz):
    for mid in np.arange(0, nm):
        for bid in np.arange(0, nb):
            
            Boost_Tab[0, bid, zid, mid] = r[idx, 0]
            Boost_Tab[1, bid, zid, mid] = r[idx, 1]
            idx = idx + 1

# Leave some redundencies in axis
Margin = 1e-5
Min = 1 - Margin
Max = 1 + Margin
zp_vec = z_vec + 1
OmB_vec[0] = OmB_vec[0] * Min
OmB_vec[-1] = OmB_vec[-1] * Max
m_vec[0] = m_vec[0] * Min
m_vec[-1] = m_vec[-1] * Max
zp_vec[0] = zp_vec[0] * Min
zp_vec[-1] = zp_vec[-1] * Max

np.savez(data_file, Boost_Tab = Boost_Tab, OmB_vec = OmB_vec, m_vec = m_vec, zp_vec = zp_vec)
print(np.shape(Boost_Tab))
print(np.max(Boost_Tab))
