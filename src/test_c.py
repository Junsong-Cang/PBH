# gcc -shared -o mdot_kernel_c.so mdot_kernel_c.c
from mdot_kernel import *
import platform, ctypes

if platform.system() == 'Darwin':
    main_path = '/Users/cangtao/cloud/GitHub/PBH/'
else:
    # running on IHEP server
    main_path = '/home/dm/watson/work/Accreting_PBH/'

c_function_lib_file = main_path + 'src/mdot_kernel_c.so'
c_function_lib = ctypes.CDLL(c_function_lib_file)

def Get_mdot_c(
    z = 10,
    m = 1000,
    OmB = 0.049,
    Use_EoR = True):
    
    xe, tk = LCDM_HyRec(z = z, Use_EoR = Use_EoR)
    # print(xe, tk)
    # Declare the function signature
    Double = ctypes.c_double
    Int = ctypes.c_int

    # specify the name of function you want
    c_function = c_function_lib.mdot_kernel_c
    
    # set input type
    c_function.argtypes = (Double, Double, Double, Double, Double, Int)
    
    # set result type
    c_function.restype = Double
    
    # Call the C function
    r = c_function(z, m, xe, tk, OmB, 1)
    
    return r

def Get_mdot_p(
    z = 10,
    m = 1000,
    OmB = 0.049,
    Use_EoR = True):
    
    r = mdot_kernel(z = z, M_PBH = m, Omega_b=OmB, Use_halo = 1, Use_EoR=Use_EoR)
    return r

reload = 1
nz = 20
nm = 20

z = np.logspace(0.5, 3, nz)
m = np.logspace(0, 4, nm)

OmB = 0.049
Use_EoR = False

if reload:
    r1 = np.zeros((nz, nm))
    r2 = np.zeros((nz, nm))
    dif = np.zeros((nz, nm))
    t1 = TimeNow()
    
    for mid in np.arange(0, nm):
        print(mid/nm)
        for zid in np.arange(0, nz):
            r1[zid, mid] = Get_mdot_c(z=z[zid], m = m[mid], OmB = OmB, Use_EoR = Use_EoR)
            # r2[zid, mid] = Get_mdot_p(z=z[zid], m = m[mid], OmB = OmB, Use_EoR = Use_EoR)
            r2[zid, mid] = 1
            dif[zid, mid] = abs(1 - r1[zid, mid]/r2[zid, mid])
    Timer(t1)
    np.savez('tmp.npz', r1 = r1, r2 = r2, dif = dif)

r = np.load('tmp.npz')
dif = r['dif']
r1 = r['r1']
print(dif.max())
print(dif.min())
# print(r1.max())
print(np.sum(dif)/400)
