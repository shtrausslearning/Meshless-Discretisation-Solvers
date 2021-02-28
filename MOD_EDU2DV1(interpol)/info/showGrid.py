# Making imports
import numpy as np
import matplotlib.pyplot as plt

# Preprocessing Input data
# physical
x,y = np.loadtxt('xyPall.dat',delimiter=',',unpack=True)
# xFF,yFF = np.loadtxt('xyPFF.dat',delimiter=',',unpack=True)
# xW,yW = np.loadtxt('xyPW.dat',delimiter=',',unpack=True)
# xBW,yBW = np.loadtxt('xyPIBW.dat',delimiter=',',unpack=True)

gcidx,gcidy = np.loadtxt('gcidxyALL.dat',delimiter=',',unpack=True)

# ghost nodes
#GCxFF,GCyFF = np.loadtxt('xyGCFF.dat',delimiter=',',unpack=True)
#GCxW,GCyW = np.loadtxt('xyGCW.dat',delimiter=',',unpack=True)
# IP nodes
#IPxFF,IPyFF = np.loadtxt('xyIPFF.dat',delimiter=',',unpack=True)
#IPxW,IPyW = np.loadtxt('xyIPW.dat',delimiter=',',unpack=True)

# specific IP's neighbours
#IPnx,IPny = np.loadtxt('IPsNEIGH.dat',delimiter=',',unpack=True)

# PHYSICAL NODES ONLY
plt.scatter(x,y,s=10,facecolors='black')

# IMAGE POINT
#plt.scatter(IPxFF,IPyFF,s=10,facecolors='purple')
#plt.scatter(IPxW,IPyW,s=10,facecolors='orange')

# GHOST CELL NODES ONLY
#plt.scatter(GCxFF,GCyFF,s=25,facecolors='none',edgecolors='purple')
#plt.scatter(GCxW,GCyW,s=25,facecolors='none',edgecolors='orange')

#plt.scatter(IPnx,IPny,s=25,facecolors='red',edgecolors='red')
plt.scatter(gcidx,gcidy,s=25,facecolors='purple',edgecolors='none')


# show near wall
plt.xlim(-0.05, 0.25)
plt.ylim((-0.045,0.1))
#plt.xlim(40, 60)
#plt.ylim(-10,10)

plt.show()