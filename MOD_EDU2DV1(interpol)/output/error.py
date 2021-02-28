import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties

font = {'family': 'sans-serif',
		'color' : 'black',
		'weight': 'bold',
		 }

label1 = 'roe'
label2 = 'rhll'
label3 = 'name3'
ncomp = 1 

# Cp OUTPUT
x,Cp  = np.loadtxt('./1/surfCp.dat', delimiter=',', unpack=True)
if( ncomp == 2 ): x,Cp2 = np.loadtxt('./2/surfCp.dat', delimiter=',', unpack=True)
if( ncomp == 3 ): x,Cp3 = np.loadtxt('./3/surfCp.dat', delimiter=',', unpack=True)
fig,ax = plt.subplots()
ax.plot(x,Cp,'o',mfc='none',mec='k',ms=4,mew=1.5,label=label1)
if( ncomp == 2 ): ax.plot(x,Cp2,'o',mfc='none',mec='b',ms=4,mew=1.5,label=label2)
if( ncomp == 3 ): ax.plot(x,Cp3,'o',mfc='none',mec='r',ms=4,mew=1.5,label=label3)
plt.xlabel(r'x coordinate')
plt.ylabel(r'-Cp')
plt.grid(True,linestyle='--')
plt.legend()
plt.savefig('./surfplots/Cp.png', bbox_inches='tight',figsize=(8,6),dpi=100)
plt.savefig('./surfplots/Cp.pdf', bbox_inches='tight',figsize=(8,6),dpi=100,format='pdf')

# CV1 OUTPUT
#x,Cp  = np.loadtxt('./1/surfcv1.dat', delimiter=',', unpack=True)
#x,Cp2 = np.loadtxt('./2/surfcv1.dat', delimiter=',', unpack=True)
#x,Cp3 = np.loadtxt('./3/surfcv1.dat', delimiter=',', unpack=True)
#fig,ax = plt.subplots()
#ax.plot(x,Cp,'o',mfc='none',mec='k',ms=4,mew=1.5,label=label1)
#ax.plot(x,Cp2,'o',mfc='none',mec='b',ms=4,mew=1.5,label=label2)
#ax.plot(x,Cp3,'o',mfc='none',mec='r',ms=4,mew=1.5,label=label3)
#plt.xlabel(r'x coordinate')
#plt.ylabel(r'conservative variable, density')
#plt.grid(True,linestyle='--')
#plt.legend()
#plt.savefig('cv1.png', bbox_inches='tight',figsize=(8,6),dpi=100)
#plt.savefig('cv1.pdf', bbox_inches='tight',figsize=(8,6),dpi=100,format='pdf')
#
## CV2 OUTPUT
#x,cv2 = np.loadtxt('surfcv2.dat', delimiter=',', unpack=True)
#fig,ax = plt.subplots()
#ax.plot(x,cv2,'o',mfc='none',mec='k',ms=4,mew=1.5,label='TLS')
#plt.xlabel(r'x coordinate')
#plt.ylabel(r'conservative variable, x momentum')
#plt.grid(True,linestyle='--')
#plt.savefig('cv2.png', bbox_inches='tight',figsize=(8,6),dpi=100)
#plt.savefig('cv2.pdf', bbox_inches='tight',figsize=(8,6),dpi=100,format='pdf')
#
## CV3 OUTPUT
#x,cv3 = np.loadtxt('surfcv3.dat', delimiter=',', unpack=True)
#fig,ax = plt.subplots()
#ax.plot(x,cv3,'o',mfc='none',mec='k',ms=4,mew=1.5,label='TLS')
#plt.xlabel(r'x coordinate')
#plt.ylabel(r'conservative variable, y momentum')
#plt.grid(True,linestyle='--')
#plt.savefig('cv3.png', bbox_inches='tight',figsize=(8,6),dpi=100)
#plt.savefig('cv3.pdf', bbox_inches='tight',figsize=(8,6),dpi=100,format='pdf')
#
## CV4 OUTPUT
#x,cv4 = np.loadtxt('surfcv4.dat', delimiter=',', unpack=True)
#fig,ax = plt.subplots()
#ax.plot(x,cv4,'o',mfc='none',mec='k',ms=4,mew=1.5,label='TLS')
#plt.xlabel(r'x coordinate')
#plt.ylabel(r'conservative variable, Total Energy')
#plt.grid(True,linestyle='--')
#plt.savefig('cv4.png', bbox_inches='tight',figsize=(8,6),dpi=100)
#plt.savefig('cv4.pdf', bbox_inches='tight',figsize=(8,6),dpi=100,format='pdf')
#
## ENTHALPY OUTPUT
#x,enth = np.loadtxt('surfenthalpy.dat', delimiter=',', unpack=True)
#fig,ax = plt.subplots()
#ax.plot(x,enth,'o',mfc='none',mec='k',ms=4,mew=1.5,label='TLS')
#plt.xlabel(r'x coordinate')
#plt.ylabel(r'Total Enthalpy')
#plt.grid(True,linestyle='--')
#plt.savefig('enthalpy.png', bbox_inches='tight',figsize=(8,6),dpi=100)
#plt.savefig('enthalpy.pdf', bbox_inches='tight',figsize=(8,6),dpi=100,format='pdf')
#
## ENTROPY OUTPUT
#x,entr = np.loadtxt('surfentropy.dat', delimiter=',', unpack=True)
#fig,ax = plt.subplots()
#ax.plot(x,entr,'o',mfc='none',mec='k',ms=4,mew=1.5,label='TLS')
#plt.xlabel(r'x coordinate')
#plt.ylabel(r'Entropy')
#plt.grid(True,linestyle='--')
#plt.savefig('entropy.png', bbox_inches='tight',figsize=(8,6),dpi=100)
#plt.savefig('entropy.pdf', bbox_inches='tight',figsize=(8,6),dpi=100,format='pdf')
#
## MACH OUTPUT
#x,mach = np.loadtxt('surfmach.dat', delimiter=',', unpack=True)
#fig,ax = plt.subplots()
#ax.plot(x,mach,'o',mfc='none',mec='k',ms=4,mew=1.5,label='TLS')
#plt.xlabel(r'x coordinate')
#plt.ylabel(r'Mach Number')
#plt.grid(True,linestyle='--')
#plt.savefig('mach.png', bbox_inches='tight',figsize=(8,6),dpi=100)
#plt.savefig('mach.pdf', bbox_inches='tight',figsize=(8,6),dpi=100,format='pdf')
