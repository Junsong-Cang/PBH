# zmax, nz = 500, 1500 should be enough
# Stats : 
# [zmax, nz, nm] : Error@z10, Error@z8
# 500, 1500, 500 : 6%, None 
# 500, 1000, 500 : 12%, 20%
# 500, 1000, 200 : 12%, 20%
# Stangely nm = 200 shows no decrease in precision

from src.main import *
from matplotlib.colors import LogNorm

t1 = TimeNow()

r1 = Get_Evo(
    fbh0 = 1e-6,
    mc = 1e1,
    sbh = 2,
    nz = 1000,
    zmax = 500,
    Fix_mdot = 0,
    Fixed_mdot = 10,
    nm = 500,
    show_status = 1)

Timer(t1)

