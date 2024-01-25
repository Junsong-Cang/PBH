from src.src_1 import *

nk = 50
Sigma = 0.5
m1 = 1
m2 = 1e3
OmB_Lost_Ratio = 0.1

LineWidth = 2
FontSize = 18

# Initialising

k1 = m2k(m2)
k2 = m2k(m1)
kc_vec = np.logspace(np.log10(k1), np.log10(k2), nk)

def fun(x):
    r = Find_PsMax(
        kc = x,
        Sigma = Sigma,
        OmB_Lost_Ratio = OmB_Lost_Ratio,
        DeltaC = 0.45,
        map_nx = 500,
        OmB0 = cosmo['OmB'],
        zmin = 8,
        zmax = 500,
        nz = 500,
        nm = 500,
        sbh_range = [-3, 5],
        show_status = True,
        Precision = 1e-2,
        Fix_mdot = 1
    )
    return r

t1 = TimeNow()
a = fun(5)
print(a)
Timer(t1)