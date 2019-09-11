import numpy as np 
import matplotlib.pylab as plt
from qcPairingModel import QCPairing


n_work = 8 #Change n_work between 4,6 and 8

measurements = np.load('measurements{}.npy'.format(n_work))


n_simulation = 4
Emax= 2
dt = 0.005
g = 1
delta=1
t= 100*dt

Pairing1 = QCPairing(n_work,n_simulation,delta=delta,g=g,dt=dt,Emax=Emax)

eigenvalues, varEigs = Pairing1.statEig(measurements,min_measure=4)
print(eigenvalues)
print(np.sqrt(varEigs))
eigdict = {}


x = measurements[:,1]
y = measurements[:,2]
idx = np.argsort(x)
x = x[idx]
y = y[idx]
for xi in x:
	eigdict[xi] = 0
for xi, yi in zip(x,y):
	eigdict[xi] += int(yi)

x = np.array(list(eigdict.keys())).astype(np.float)
idx = np.argsort(x)
y = np.array(list(eigdict.values())).astype(np.int)
x = x[idx]
y = y[idx]

plt.plot(x,y)

plt.xlabel('Eigenvalue')
plt.ylabel('Times measured')
plt.show()


