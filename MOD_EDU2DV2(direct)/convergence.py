import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties

iter,L2R = np.loadtxt('convergence.dat', delimiter=',', unpack=True)

fig,ax = plt.subplots()
plt.setp(ax.spines.values(), linewidth=3)

ax.plot(iter,L2R ,linewidth=2, color='black')

plt.xlabel(r'iteration')
plt.ylabel(r'L2N density residual')
plt.grid(True,linestyle='--')

#plt.savefig('outD.png', bbox_inches='tight',dpi=300,quality=100)
plt.show()
