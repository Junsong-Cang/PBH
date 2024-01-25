from src.main import *

reload = 1
m = 1e-3
OmB = 0.005

nz = 50
z = np.linspace(8, 39, nz)
EoR = 1

# ----
if reload:
    t1 = TimeNow()
    r1 = np.linspace(0, 1, nz)
    r2 = np.linspace(0, 1, nz)
    
    for idx in np.arange(0, nz):
        r1[idx] = BoostFactor(z = z[idx], m = m, OmB = OmB, Use_EoR = EoR, ncpu = 12, Use_Interp = 0)
        r2[idx] = BoostFactor(z = z[idx], m = m, OmB = OmB, Use_EoR = EoR, ncpu = 12, Use_Interp = 1)
        print(idx/nz)
    np.savez('tmp.npz', r1 = r1, r2 = r2)
    Timer(t1)

r = np.load('tmp.npz')
r1 = r['r1']
r2 = r['r2']

plt.semilogy(z, r1, 'k')
plt.semilogy(z, r2, 'r')
plt.show()
