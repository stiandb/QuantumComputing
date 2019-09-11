from hamiltonian import *
import numpy as np
H = hamiltonian(1,2,1,1)
eigvals,eigvecs = np.linalg.eigh(H[0])

print('Eigenvalues for 1 pair', eigvals)

H = hamiltonian(2,2,1,1)
eigvals,eigvecs = np.linalg.eigh(H[0])
print('Eigenvalues for 2 pairs', eigvals)
