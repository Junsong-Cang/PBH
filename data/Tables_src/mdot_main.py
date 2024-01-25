nm = 300
nz = 150
z1 = 6.99
z2 = 1000.1
m1 = 1e-3
m2 = 1e5

reload = 0
ncpu = 12
EoR_z = 36

swap_file = 'swap/mdot_tab_main_swap.npz'
data_file = '/Users/cangtao/cloud/GitHub/PBH/data/mdot_tab.npz'

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

def model(id = 1):

    z = params[0, id]
    m = params[1, id]
    r = []
    r1 = mdot_kernel(z = z, M_PBH =m, Use_EoR = False)
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

r = np.load(swap_file)['r']
mdot_vec = np.zeros((2, nz, nm))

idx = 0
for zid in np.arange(0,nz):
    for mid in np.arange(0,nm):
        mdot_vec[0, zid, mid] = r[idx, 0]
        mdot_vec[1, zid, mid] = r[idx, 1]
        idx_ = round(r[idx, 2])
        if not (idx == idx_):
            raise Exception
        idx = idx + 1

# Saving data
np.savez(data_file, mdot_vec = mdot_vec, z = z_vec, m = m_vec)

# Try print things for c
c_file = '/Users/cangtao/cloud/GitHub/PBH/src/c/mdot_tab.h'

#nz = 10
#nm = 5
#x = 1.3e3*np.ones((nz, nm))

F=open(c_file,'w')

if np.min(mdot_vec) < 1e-200:
    print('#define ok_to_use_log_mdot 0', file = F)
else:
    print('#define ok_to_use_log_mdot 1', file = F)

print('#define nz ' + str(int(nz)), file = F)
print('#define nm ' + str(int(nm)), file = F)
print('#define z_min ' + "{0:.10E}".format(z_vec[0]), file = F)
print('#define z_max ' + "{0:.10E}".format(z_vec[-1]), file = F)
print('#define m_min ' + "{0:.10E}".format(m_vec[0]), file = F)
print('#define m_max ' + "{0:.10E}".format(m_vec[-1]), file = F)
print('', file = F)

# Print z_axis

s = 'double z_axis[nz] = {'
for idx in np.arange(0,nz):
    if idx == nz-1:
        s = s + "{0:.8E}".format(z_vec[idx]) + '};'
    else:
        s = s + "{0:.8E}".format(z_vec[idx]) + ', '
print(s, file = F)

s = 'double m_axis[nm] = {'
for idx in np.arange(0, nm):
    if idx == nm-1:
        s = s + "{0:.8E}".format(m_vec[idx]) + '};'
    else:
        s = s + "{0:.8E}".format(m_vec[idx]) + ', '
print(s, file = F)
print('', file = F)


# Print LCDM
x = mdot_vec[0,:,:]
print('double mdot_vec_LCDM[nz][nm] = {', file = F)
for zid in np.arange(0,nz):
    s = '   {'
    for mid in np.arange(0,nm):
        x_ = x[zid, mid]
        if mid == nm-1:
            s_ = "{0:.8E}".format(x_)
        else:
            s_ = "{0:.8E}".format(x_) + ', '
        s = s + s_
    if zid == nz-1:
        s = s + '}'
    else:
        s = s + '},'
    print(s, file = F)
print('};', file = F)
print('', file = F)

# Print EoR
x = mdot_vec[1,:,:]
print('double mdot_vec_EoR['+str(int(nz))+']['+str(int(nm))+'] = {', file = F)
for zid in np.arange(0,nz):
    s = '   {'
    for mid in np.arange(0,nm):
        x_ = x[zid, mid]
        if mid == nm-1:
            s_ = "{0:.8E}".format(x_)
        else:
            s_ = "{0:.8E}".format(x_) + ', '
        s = s + s_
    if zid == nz-1:
        s = s + '}'
    else:
        s = s + '},'
    print(s, file = F)
print('};', file = F)
print('', file = F)

F.close()
