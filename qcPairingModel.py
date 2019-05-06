import qiskit as qk
import matplotlib.pylab as plt
import numpy as np
qk.IBMQ.load_accounts()


def H0(n_simulation,control_qubit,delta,dt):
	for p in range(2,n_simulation):
		qz.crz(2*dt*(p-1),qb[control_qubit],qb[p-1+n_work])

def H1(n_simulation,control_qubit,g,dt):
	for p in range(1,n_simulation,2):
		for q in range(p,n_simulation,2):
			if p == q:
				theta = 2*(1/8)*g*dt
				qz.crz(theta,qb[control_qubit],qb[p-1+n_work])
				qz.crz(theta,qb[control_qubit],qb[p+n_work])
				qz.cx(qb[p],qb[n_qubits-1])
				qz.cx(qb[p-1],qb[n_qubits-1])
				qz.crz(theta,qb[control_qubit],qb[n_qubits-1])
				qz.cx(qb[p],qb[n_qubits-1])
				qz.cx(qb[p-1],qb[n_qubits-1])
			else:
				theta = 2*(1/16)*g*dt
				#FIRST TERM:
				qz.h(qb[p-1+n_work])
				qz.h(qb[p+n_work])
				qz.h(qb[q-1+n_work])
				qz.h(qb[q+n_work])
				qz.cx(qb[p-1+n_work],qb[n_qubits-1])
				qz.cx(qb[p+n_work],qb[n_qubits-1])
				qz.cx(qb[q-1+n_work],qb[n_qubits-1])
				qz.cx(qb[q+n_work],qb[n_qubits-1])
				qz.crz(theta,qb[control_qubit],qb[n_qubits-1])
				qz.cx(qb[p-1+n_work],qb[n_qubits-1])
				qz.cx(qb[p+n_work],qb[n_qubits-1])
				qz.cx(qb[q-1+n_work],qb[n_qubits-1])
				qz.cx(qb[q+n_work],qb[n_qubits-1])
				qz.h(qb[p-1+n_work])
				qz.h(qb[p+n_work])
				qz.h(qb[q-1+n_work])
				qz.h(qb[q+n_work])
				############
				#SECOND TERM:
				qz.h(qb[p-1+n_work])
				qz.h(qb[p+n_work])
				qz.rz(np.pi/2,qb[q-1+n_work])
				qz.h(qb[q-1+n_work])
				qz.rz(np.pi/2,qb[q+n_work])
				qz.h(qb[q+n_work])
				qz.cx(qb[p-1+n_work],qb[n_qubits-1])
				qz.cx(qb[p+n_work],qb[n_qubits-1])
				qz.cx(qb[q-1+n_work],qb[n_qubits-1])
				qz.cx(qb[q+n_work],qb[n_qubits-1])
				qz.crz(-theta,qb[control_qubit],qb[n_qubits-1])
				qz.cx(qb[p-1+n_work],qb[n_qubits-1])
				qz.cx(qb[p+n_work],qb[n_qubits-1])
				qz.cx(qb[q-1+n_work],qb[n_qubits-1])
				qz.cx(qb[q+n_work],qb[n_qubits-1])
				qz.h(qb[p-1+n_work])
				qz.h(qb[p+n_work])
				qz.h(qb[q-1+n_work])
				qz.rz(-np.pi/2,qb[q-1+n_work])
				qz.h(qb[q+n_work])
				qz.rz(-np.pi/2,qb[q+n_work])
				###########
				#THIRD TERM:
				qz.h(qb[p-1+n_work])
				qz.rz(np.pi/2,qb[p+n_work])
				qz.h(qb[p+n_work])
				qz.h(qb[q-1+n_work])
				qz.rz(np.pi/2,qb[q+n_work])
				qz.h(qb[q+n_work])
				qz.cx(qb[p-1+n_work],qb[n_qubits-1])
				qz.cx(qb[p+n_work],qb[n_qubits-1])
				qz.cx(qb[q-1+n_work],qb[n_qubits-1])
				qz.cx(qb[q+n_work],qb[n_qubits-1])
				qz.crz(theta,qb[control_qubit],qb[n_qubits-1])
				qz.cx(qb[p-1+n_work],qb[n_qubits-1])
				qz.cx(qb[p+n_work],qb[n_qubits-1])
				qz.cx(qb[q-1+n_work],qb[n_qubits-1])
				qz.cx(qb[q+n_work],qb[n_qubits-1])
				qz.h(qb[p-1+n_work])
				qz.h(qb[p+n_work])
				qz.rz(-np.pi/2,qb[p+n_work])
				qz.h(qb[q-1+n_work])
				qz.h(qb[q+n_work])
				qz.rz(-np.pi/2,qb[q+n_work])
				###########
				#FOURTH TERM
				qz.h(qb[p-1+n_work])
				qz.rz(np.pi/2,qb[p+n_work])
				qz.h(qb[p+n_work])
				qz.rz(np.pi/2,qb[q-1+n_work])
				qz.h(qb[q-1+n_work])
				qz.h(qb[q+n_work])
				qz.cx(qb[p-1+n_work],qb[n_qubits-1])
				qz.cx(qb[p+n_work],qb[n_qubits-1])
				qz.cx(qb[q-1+n_work],qb[n_qubits-1])
				qz.cx(qb[q+n_work],qb[n_qubits-1])
				qz.crz(theta,qb[control_qubit],qb[n_qubits-1])
				qz.cx(qb[p-1+n_work],qb[n_qubits-1])
				qz.cx(qb[p+n_work],qb[n_qubits-1])
				qz.cx(qb[q-1+n_work],qb[n_qubits-1])
				qz.cx(qb[q+n_work],qb[n_qubits-1])
				qz.h(qb[p-1+n_work])
				qz.h(qb[p+n_work])
				qz.rz(-np.pi/2,qb[p+n_work])
				qz.h(qb[q-1+n_work])
				qz.rz(-np.pi/2,qb[q-1+n_work])
				qz.h(qb[q+n_work])
				###########
				#FIFTH TERM
				qz.rz(np.pi/2,qb[p-1+n_work])
				qz.h(qb[p-1+n_work])
				qz.h(qb[p+n_work])
				qz.h(qb[q-1+n_work])
				qz.rz(np.pi/2,qb[q+n_work])
				qz.h(qb[q+n_work])
				qz.cx(qb[p-1+n_work],qb[n_qubits-1])
				qz.cx(qb[p+n_work],qb[n_qubits-1])
				qz.cx(qb[q-1+n_work],qb[n_qubits-1])
				qz.cx(qb[q+n_work],qb[n_qubits-1])
				qz.crz(theta,qb[control_qubit],qb[n_qubits-1])
				qz.cx(qb[p-1+n_work],qb[n_qubits-1])
				qz.cx(qb[p+n_work],qb[n_qubits-1])
				qz.cx(qb[q-1+n_work],qb[n_qubits-1])
				qz.cx(qb[q+n_work],qb[n_qubits-1])
				qz.h(qb[p-1+n_work])
				qz.rz(-np.pi/2,qb[p-1+n_work])
				qz.h(qb[p+n_work])
				qz.h(qb[q-1+n_work])
				qz.h(qb[q+n_work])
				qz.rz(-np.pi/2,qb[q+n_work])
				##########
				#SIXTH TERM:
				qz.rz(np.pi/2,qb[p-1+n_work])
				qz.h(qb[p-1+n_work])
				qz.h(qb[p+n_work])
				qz.rz(np.pi/2,qb[q-1+n_work])
				qz.h(qb[q-1+n_work])
				qz.h(qb[q+n_work])
				qz.cx(qb[p-1+n_work],qb[n_qubits-1])
				qz.cx(qb[p+n_work],qb[n_qubits-1])
				qz.cx(qb[q-1+n_work],qb[n_qubits-1])
				qz.cx(qb[q+n_work],qb[n_qubits-1])
				qz.crz(theta,qb[control_qubit],qb[n_qubits-1])
				qz.cx(qb[p-1+n_work],qb[n_qubits-1])
				qz.cx(qb[p+n_work],qb[n_qubits-1])
				qz.cx(qb[q-1+n_work],qb[n_qubits-1])
				qz.cx(qb[q+n_work],qb[n_qubits-1])
				qz.h(qb[p-1+n_work])
				qz.rz(-np.pi/2,qb[p-1+n_work])
				qz.h(qb[p+n_work])
				qz.h(qb[q-1+n_work])
				qz.rz(-np.pi/2,qb[q-1+n_work])
				qz.h(qb[q+n_work])
				#######################
				#SEVENTH TERM
				qz.rz(np.pi/2,qb[p-1+n_work])
				qz.h(qb[p-1+n_work])
				qz.rz(np.pi/2,qb[p+n_work])
				qz.h(qb[p+n_work])
				
				qz.h(qb[q-1+n_work])
				
				qz.h(qb[q+n_work])
				qz.cx(qb[p-1+n_work],qb[n_qubits-1])
				qz.cx(qb[p+n_work],qb[n_qubits-1])
				qz.cx(qb[q-1+n_work],qb[n_qubits-1])
				qz.cx(qb[q+n_work],qb[n_qubits-1])
				qz.crz(-theta,qb[control_qubit],qb[n_qubits-1])
				qz.cx(qb[p-1+n_work],qb[n_qubits-1])
				qz.cx(qb[p+n_work],qb[n_qubits-1])
				qz.cx(qb[q-1+n_work],qb[n_qubits-1])
				qz.cx(qb[q+n_work],qb[n_qubits-1])
				qz.h(qb[p-1+n_work])
				qz.rz(-np.pi/2,qb[p-1+n_work])
				qz.h(qb[p+n_work])
				qz.rz(-np.pi/2,qb[p+n_work])
				qz.h(qb[q-1+n_work])
				
				qz.h(qb[q+n_work])
				##############
				#EIGTH TERM:
				qz.rz(np.pi/2,qb[p-1+n_work])
				qz.h(qb[p-1+n_work])
				qz.rz(np.pi/2,qb[p+n_work])
				qz.h(qb[p+n_work])
				qz.rz(np.pi/2,qb[q-1+n_work])
				qz.h(qb[q-1+n_work])
				qz.rz(np.pi/2,qb[q+n_work])
				qz.h(qb[q+n_work])
				qz.cx(qb[p-1+n_work],qb[n_qubits-1])
				qz.cx(qb[p+n_work],qb[n_qubits-1])
				qz.cx(qb[q-1+n_work],qb[n_qubits-1])
				qz.cx(qb[q+n_work],qb[n_qubits-1])
				qz.crz(theta,qb[control_qubit],qb[n_qubits-1])
				qz.cx(qb[p-1+n_work],qb[n_qubits-1])
				qz.cx(qb[p+n_work],qb[n_qubits-1])
				qz.cx(qb[q-1+n_work],qb[n_qubits-1])
				qz.cx(qb[q+n_work],qb[n_qubits-1])
				qz.h(qb[p-1+n_work])
				qz.rz(-np.pi/2,qb[p-1+n_work])
				qz.h(qb[p+n_work])
				qz.rz(-np.pi/2,qb[p+n_work])
				qz.h(qb[q-1+n_work])
				qz.rz(-np.pi/2,qb[q-1+n_work])
				qz.h(qb[q+n_work])
				qz.rz(-np.pi/2,qb[q+n_work])



n_work = 8
n_simulation = 4
n_qubits = n_work + n_simulation + 1


qb = qk.QuantumRegister(n_qubits)
c = qk.ClassicalRegister(n_qubits)
qz = qk.QuantumCircuit(qb,c)
for i in range(n_work):
	qz.h(qb[i])

for i in range(n_simulation):
	qz.h(qb[i+n_work])

H0(n_simulation,0,1,0.1)
H1(n_simulation,0,1,0.1)

qz.measure(qb,c)

job = qk.execute(qz, backend = qk.Aer.get_backend('qasm_simulator'), shots=10000)
result = job.result()

# Print the result
print(result.get_counts(qz))
