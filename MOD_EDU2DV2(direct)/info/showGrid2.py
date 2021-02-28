# Making imports
import numpy as np
import matplotlib.pyplot as plt

# Preprocessing Input data

# physical
#x,y = np.loadtxt('xyPall.dat',delimiter=',',unpack=True)
# xFF,yFF = np.loadtxt('xyPFF.dat',delimiter=',',unpack=True)
# xW,yW = np.loadtxt('xyPW.dat',delimiter=',',unpack=True)
# xBW,yBW = np.loadtxt('xyPIBW.dat',delimiter=',',unpack=True)

# GCID new
#gcidx,gcidy = np.loadtxt('gcidnew1.dat',delimiter=',',unpack=True)
#gcidx2,gcidy2 = np.loadtxt('gcidnew2.dat',delimiter=',',unpack=True)

#specific node neighbours
#xnei,ynei = np.loadtxt('wallneighb.dat',delimiter=',',unpack=True)

#output Aij Coefficients
aij,aji = np.loadtxt('lsqCoeffAij.dat',delimiter=',',unpack=True)



# PHYSICAL NODES ONLY X,Y
#plt.scatter(x,y,s=10,facecolors='black')
#plt.scatter(gcidx,gcidy,s=10,facecolors='blue')
#plt.scatter(gcidx2,gcidy2,s=10,facecolors='red')
#plt.scatter(xnei,ynei,s=10,facecolors='orange')

plt.scatter(aij,aji)

# show near wall
plt.xlim(-0.05, 0.25)
plt.ylim((-0.045,0.1))
#plt.xlim(0.9,1.1)
#plt.ylim(-0.1,0.1)
#plt.xlim(0.99,1.01)
#plt.ylim(-0.01,0.01)
#plt.xlim(40, 60)
#plt.ylim(-10,10)

plt.show()