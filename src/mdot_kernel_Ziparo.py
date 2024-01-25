# Zip's mdot module
from PyLab import Read_Curve
import numpy as np
import matplotlib.pyplot as plt
# import costanti as ct

DataPath = '/Users/cangtao/cloud/Library/PyLab/Curve_Data/'

def vrel2(z):
	zp = 1+z
	f = DataPath + '0709.0524.Fig2.vrel.txt'
	zp_axis, vrel_axis = Read_Curve(
		File = f,
		nx = 1000,
		model = 2,
		Convert_x = 1,
		Convert_y = 1
		)
	vrel_axis = vrel_axis
	vr = np.interp(x = zp, xp = zp_axis, fp = vrel_axis)
	return vr

#------------------ Effective velocity ------------------#
def vrel2_zip(z):
	if (z<200):
		a = ((1+z)/10)**(1/2.4)*2
		return a
	else:
		a=(200/(1+z))**(1/2.4)*((1+200)/10)**(1/2.4)*2
		return a

def cs(z):
	b=1.72
	z_dec=130
	a = 5.74*((1+z)/1000)**(1./2.)*(((1+z_dec)/(1+z))**b+1)**(-(1/(2*b)))
	return a

def mach(z):
	mach=vrel2(z)/cs(z)
	return mach

def veff(z):
	b = mach(z)
	if b>1:
		a = cs(z)*(16/np.sqrt(2*np.pi)*mach(z)**3)**(1/6)
		return a
	else:
		a = cs(z)*(1+mach(z)**2)**(1/2)
		return a
#------------------ Effective velocity ------------------#

#------------------ Accretion ratio ------------------#
def beta(z, mbh):
	xe=1e-3
	return (mbh/1e+4)*((1+z)/1000.)**(3/2)*(veff(z)/5.74)**(-3)*(0.257+1.45*(xe/0.01)*((1+z)/1000)**(5./2.))

def xc(z, mbh):
	a = (-1+(1+beta(z,mbh))**(1./2.))/beta(z,mbh)
	return a

def l(z, mbh):
	return xc(z, mbh)**2*np.exp(4.5/(3+beta(z, mbh)**(0.75))) 

def ratio(z, mbh):
	a = (0.0018)*l(z, mbh)*((1+z)/1000)**3*mbh*(veff(z)/5.74)**(-3)
	return a

#------------------ Accretion ratio ------------------#

zmin=3
zmax=4e+3
z = np.linspace(zmin, zmax, 1000)
mass=[1, 10, 1e+3, 1e+4, 1e+5]
Lumi= np.zeros_like(z)
Lumi2= np.zeros_like(z)
Lumi3= np.zeros_like(z)

# Paper results
f = DataPath + '0709.0524.Fig3a.top1.txt'
zp, r2 = Read_Curve(
        File = f,
        nx = 100,
        model = 2,
        Convert_x = 1,
        Convert_y = 1
        )

for j in range (0, len(mass)):
	PBHMass=mass[j]
	for i in range (0, len(z)):
		Lumi[i]=ratio(z[i], PBHMass)
	if j == 4:
		plt.plot(z, Lumi, 'k')
plt.plot(zp - 1, r2, 'r')
plt.xlim(zmin, zmax)
#plt.ylim(5*1e-7, 10)
plt.xscale('log')
plt.yscale('log')
# plt.plot(z, Lumi)
# plt.plot(z, Lumi2, 'k')

print(Lumi2)

plt.show()

	

