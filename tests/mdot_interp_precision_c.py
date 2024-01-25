from src.main import *

reload = 0
nz = 1051
nm = 541
Use_EoR = 1
Use_Edd_Limit = 0

z_vec = np.logspace(1, 3, nz)
m_vec = np.logspace(-3, 4, nm)

r0 = np.zeros((nz, nm))
r1 = np.zeros((nz, nm))

if reload:
    t1 = TimeNow()
    for zid in np.arange(0, nz):
        for mid in np.arange(0, nm):
            r0[zid, mid] = Get_mdot(z = z_vec[zid], m = m_vec[mid], Use_EoR = Use_EoR, Use_Edd_Limit = Use_Edd_Limit, Use_C = 0)
    Timer(t1)
    t1 = TimeNow()
    for zid in np.arange(0, nz):
        for mid in np.arange(0, nm):
            r1[zid, mid] = Get_mdot(z = z_vec[zid], m = m_vec[mid], Use_EoR = Use_EoR, Use_Edd_Limit = Use_Edd_Limit, Use_C = 1)
    Timer(t1)
    np.savez('tmp.npz', r0 = r0, r1 = r1)

r = np.load('tmp.npz')
r0 = r['r0']
r1 = r['r1']

dif = np.zeros((nz, nm))

for zid in np.arange(0, nz):
    for mid in np.arange(0, nm):
        r0_ = r0[zid,mid]
        r1_ = r1[zid,mid]
        dif[zid, mid] = np.abs(r0_ - r1_)/(min(r0_, r1_))

dif_avg = np.sum(dif)/(nz*nm)

print('dif_max = ', dif.max())
print('dif_min = ', dif.min())
print('dif_avg = ', dif_avg)

idx = np.argmax(dif)
zid = int(np.floor(idx/nm))
mid = idx - nm * zid
xx = dif[zid, mid]
print('recovered dif = ', xx)
m_ = m_vec[mid]
z_ = z_vec[zid]
print('bogey : z = ', z_, ', m = ', m_)

d0 = Get_mdot(z = z_vec[zid], m = m_vec[mid], Use_EoR = Use_EoR, Use_Edd_Limit = Use_Edd_Limit, Do_not_Interpolate = 1, Use_C = 0)
d1 = Get_mdot(z = z_vec[zid], m = m_vec[mid], Use_EoR = Use_EoR, Use_Edd_Limit = Use_Edd_Limit, Do_not_Interpolate = 0, Use_C = 0)
d2 = Get_mdot(z = z_vec[zid], m = m_vec[mid], Use_EoR = Use_EoR, Use_Edd_Limit = Use_Edd_Limit, Do_not_Interpolate = 0, Use_C = 1)

print(d0, d1, d2)
