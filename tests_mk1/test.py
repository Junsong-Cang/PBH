# test that all curves are corerctly extracted

from PyLab import *

DataPath = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/'

f = DataPath + '0709.0524.Fig2.cs.txt'
f = DataPath + '0709.0524.Fig2.vrel.txt'
f = DataPath + '0709.0524.Fig2.veff_A.txt'
f = DataPath + '0709.0524.Fig2.veff_B.txt'

f = DataPath + '0709.0524.Fig3a.top1.txt'
f = DataPath + '0709.0524.Fig3a.top2.txt'
f = DataPath + '0709.0524.Fig3a.top3.txt'
f = DataPath + '0709.0524.Fig3a.top4.txt'

f = DataPath + '0709.0524.Fig3b.top1.txt'
f = DataPath + '0709.0524.Fig3b.top2.txt'
f = DataPath + '0709.0524.Fig3b.top3.txt'
f = DataPath + '0709.0524.Fig3b.top4.txt'

f = DataPath + '0709.0524.Fig4a.top1.txt'
f = DataPath + '0709.0524.Fig4a.top2.txt'
f = DataPath + '0709.0524.Fig4a.top3.txt'
f = DataPath + '0709.0524.Fig4a.top4.txt'
f = DataPath + '0709.0524.Fig4a.dotted.txt'

f = DataPath + '0709.0524.Fig4b.top1.txt'
f = DataPath + '0709.0524.Fig4b.top2.txt'
f = DataPath + '0709.0524.Fig4b.top3.txt'
f = DataPath + '0709.0524.Fig4b.top4.txt'
f = DataPath + '0709.0524.Fig4b.dotted.txt'

f = DataPath + '2003.02778.fig1.m1e5.txt'
f = DataPath + '2003.02778.fig1.m1e4.txt'
f = DataPath + '2003.02778.fig1.m3e2.txt'
f = DataPath + '2003.02778.fig1.m1e2.txt'
f = DataPath + '2003.02778.fig1.m10.txt'
f = DataPath + '2003.02778.fig1.m1.txt'

x, y = Read_Curve(
    File = f,
    model = 2,
    Convert_x = 1,
    Convert_y = 1
)

plt.loglog(x, y)
plt.show()
