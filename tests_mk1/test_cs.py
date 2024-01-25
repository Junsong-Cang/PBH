# test that all curves are corerctly extracted

from PyLab import *

DataPath = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/'

f = DataPath + '0709.0524.Fig2.cs.txt'

zp, cs1 = Read_Curve(
    File = f,
    model = 2,
    Convert_x = 1,
    Convert_y = 1
)

cs1 = cs1*1e3

xe, Tk = LCDM_HyRec(z = zp-1)
cs2 = 5700 * (Tk/2730)**0.5

plt.loglog(zp, cs1, 'k')
plt.loglog(zp, cs2, 'r')
plt.loglog(zp, cs3, 'b')

plt.show()
