from src.main import *

LineWidth = 2
FontSize = 15

reload = 0

nz = 29
nm = 23
m1 = 1e2
m2 = 1e5
z1 = 7.1
z2 = 999.9
Use_EoR = 1
ncpu = 12
datafile = 'data/mdot_interp_2D.npz'

# ---- starting ----
# Stats:
# m = [1e-3, 1e0] : avg_dif = 0.003, max_dif = 0.06, max_dif_loc: z = 12.73275, m = 1

diff = lambda x1, x2: np.abs(x2 - x1)/(min(x1, x2))

from matplotlib.colors import LogNorm
z = np.logspace(np.log10(z1), np.log10(z2), nz)
m = np.logspace(np.log10(m1), np.log10(m2), nm)
zids = np.arange(0, nz)

def Get_dif(id):
    z_ = z[id]
    r = np.linspace(0, 1, nm)
    status = id/nz
    for mid in np.arange(0, nm):
        m_ = m[mid]
        x1 = Get_mdot(z = z_, m = m_, Use_Edd_Limit = 0, Use_EoR = Use_EoR)
        x2 = Get_mdot(z = z_, m = m_, Use_Edd_Limit = 0, Do_not_Interpolate=True, Use_EoR = Use_EoR)
        dif = diff(x1, x2)
        
        print('dif = ', dif)
        r[mid] = dif
    print('--------status--------', status)
    return r

if reload:
    t1 = TimeNow()
    dif = Parallel(n_jobs=ncpu)(delayed(Get_dif)(x) for x in zids)
    Timer(t1)
    np.savez(datafile, dif = dif)

dif = np.load(datafile)['dif']
dif_min = dif.min()
dif_max = dif.max()
avg = np.sum(dif)/nz/nm

# Print preliminary stats
print('avg_dif = ', avg)
print('dif_min = ', dif.min())
print('dif_max = ', dif_max)

# dif index : [zid, mid]

idx = np.argmax(dif)

zid = int(np.floor(idx/nm))
mid = idx - nm * zid
print('max_dif_index = ', idx)
print('zid = ', zid)
print('mid = ', mid)

#zid = 2
#mid = 22

xx = dif[zid, mid]

print('recovered dif = ', xx)

m_ = m[mid]
z_ = z[zid]

print('bogey coordinates [z , m] = [', z_, ', ', m_, ']')

x1 = Get_mdot(z = z_, m = m_, Use_EoR = Use_EoR, Use_Edd_Limit = 0)
x2 = Get_mdot(z = z_, m = m_, Do_not_Interpolate=True, Use_EoR = Use_EoR, Use_Edd_Limit = 0)
d = diff(x1, x2)
print('Re-computed dif : ')
print('interp = ', x1, ', no_interp = ', x2, ', dif =', d)

# Load table
Tab = Interpolation_Table['mdot_data']
m_axis = Tab['m']
if Use_EoR:
    idx = 1
else:
    idx = 0
mdot_tab = Tab['mdot_vec'][idx,:,:]
zp_axis = Tab['z'] + 1

zid1 = Find_Index(x = z_ + 1, x_axis = zp_axis)
mid1 = Find_Index(x = m_, x_axis = m_axis)
zid2 = zid1 + 1
mid2 = mid1 + 1

z1 = zp_axis[zid1] - 1
z2 = zp_axis[zid2] - 1
m1 = m_axis[mid1]
m2 = m_axis[mid2]

f11 = mdot_tab[zid1, mid1]
f12 = mdot_tab[zid1, mid2]
f21 = mdot_tab[zid2, mid1]
f22 = mdot_tab[zid2, mid2]
print('--Interp Table Coordinates--')
print('[z1, m1] = ', [z1, m1], ', f11 = ', f11)
print('[z1, m1] = ', [z1, m2], ', f12 = ', f12)
print('[z2, m1] = ', [z2, m1], ', f21 = ', f21)
print('[z2, m2] = ', [z2, m2], ', f22 = ', f22)

print('z_axis range:', [zp_axis[0]-1, zp_axis[-1] - 1])

print('--stats-- @ [z1, m1]:')
x1 = Get_mdot(z = z1, m = m1, Use_EoR = Use_EoR, Use_Edd_Limit = 0)
x2 = Get_mdot(z = z1, m = m1, Do_not_Interpolate=True, Use_EoR = Use_EoR, Use_Edd_Limit = 0)
d = diff(x1, x2)
print('interp = ', x1, ', no_interp = ', x2, ', dif =', d)

print('--stats-- @ [z1, m2]:')
x1 = Get_mdot(z = z1, m = m2, Use_EoR = 1, Use_Edd_Limit = 0)
x2 = Get_mdot(z = z1, m = m2, Do_not_Interpolate=True, Use_EoR = Use_EoR, Use_Edd_Limit = 0)
d = diff(x1, x2)
print('interp = ', x1, ', no_interp = ', x2, ', dif =', d)

print('--stats-- @ [z2, m1]:')

x1 = Get_mdot(z = z2, m = m1, Use_EoR = 1, Use_Edd_Limit = 0)
x2 = Get_mdot(z = z2, m = m1, Do_not_Interpolate=True, Use_EoR = Use_EoR, Use_Edd_Limit = 0)

d = diff(x1, x2)

print('interp = ', x1, ', no_interp = ', x2, ', dif =', d)

print('--stats-- @ [z2, m2]:')
x1 = Get_mdot(z = z2, m = m2, Use_EoR = 1, Use_Edd_Limit = 0)
x2 = Get_mdot(z = z2, m = m2, Do_not_Interpolate=True, Use_EoR = Use_EoR, Use_Edd_Limit = 0)
d = diff(x1, x2)
print('interp = ', x1, ', no_interp = ', x2, ', dif =', d)

# ----Get plot----


plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
m = np.arange(0,len(m))
z = np.arange(0,len(z))

m, z = np.meshgrid(m, z)

fig,ax = plt.subplots()
PlotData = 1+dif
c=ax.pcolor(m, z, PlotData, cmap='jet',norm = LogNorm(vmin=PlotData.min(), vmax=PlotData.max()))
# c=ax.pcolor(m, z, PlotData, cmap='jet')

plt.xlabel('$m$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$z$',fontsize=FontSize,fontname='Times New Roman')
#plt.xscale('log')
#plt.yscale('log')

plt.xticks(size=FontSize)
plt.yticks(size=FontSize)

clb = plt.colorbar(c)
clb.ax.set_title('dif')

plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.png',dpi=200)

