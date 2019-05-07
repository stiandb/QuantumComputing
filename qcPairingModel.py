import qiskit as qk
import matplotlib.pylab as plt
import numpy as np
qk.IBMQ.load_accounts()
E_max = 1e16
def H0(n_simulation,control_qubit,delta,dt):
	for p in range(2,n_simulation):
		qz.crz(2*dt*(p-1),qb[control_qubit],qb[p-1+n_work])

	qz.cu1(E_max*dt,qb[control_qubit],qb[n_work])
	qz.x(qb[n_work])
	qz.cu1(E_max*dt,qb[control_qubit],qb[n_work])
	qz.x(qb[n_work])

	qz.cu1(-dt*delta*0.5*(n_simulation-1)*n_simulation,qb[control_qubit],qb[n_work])
	qz.x(qb[n_work])
	qz.cu1(-dt*delta*0.5*(n_simulation-1)*n_simulation,qb[control_qubit],qb[n_work])
	qz.x(qb[n_work])

def H1(n_simulation,control_qubit,g,dt):
	for p in range(1,n_simulation,2):
		for q in range(p,n_simulation,2):
			if p == q:
				theta = 2*(1/8)*g*dt
				qz.cu1(-theta/2,qb[control_qubit],qb[n_work])
				qz.x(qb[n_work])
				qz.cu1(-theta/2,qb[control_qubit],qb[n_work])
				qz.x(qb[n_work])

				qz.crz(theta,qb[control_qubit],qb[p-1+n_work])
				qz.crz(theta,qb[control_qubit],qb[p+n_work])
				
				qz.cx(qb[p],qb[n_qubits-1])
				qz.cx(qb[p-1],qb[n_qubits-1])
				qz.crz(-theta,qb[control_qubit],qb[n_qubits-1])
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


def inverse_fourier():
	for qc in range(int(n_work/2)):
		qz.swap(qb[qc],qb[n_work-qc-1])
	
	for cq in range(n_work):
		for i in range(cq):
			qz.cu1(-2*np.pi/(2**(1+cq-i)),qb[i],qb[cq])
		qz.h(qb[cq])

n_work = 4
n_simulation = 2
n_qubits = n_work + n_simulation + 1


qb = qk.QuantumRegister(n_qubits)
c = qk.ClassicalRegister(n_qubits)
qz = qk.QuantumCircuit(qb,c)
delta = 1
g = 0.001
dt = 2*np.pi/(E_max + 100)
t = dt
steps = int(t/dt)

for i in range(n_simulation):
	qz.x(qb[n_work+i])

for cq in range(n_work):
	qz.h(qb[cq])
for cq in range(n_work):
	for j in range(steps):
		H0(n_simulation,cq,delta,(2**cq)*dt)
		H1(n_simulation,cq,g,(2**cq)*dt)

inverse_fourier()

qz.measure(qb,c)

shots = 1000
job = qk.execute(qz, backend = qk.Aer.get_backend('qasm_simulator'), shots=shots)
result = job.result()

result = result.get_counts(qz)
res_decimal = {}
res_energy = {}
for key,value in result.items():
	key = key[0:n_work]
	decimal = 0
	for i,bit in enumerate(key):
		decimal += int(bit)*2**(-i-1)
	if value != 0:
		res_decimal[str(decimal)] = 0
		res_energy[str(-decimal*2*np.pi/dt)] = 0
for key,value in result.items():
	key = key[0:n_work]
	decimal = 0
	for i,bit in enumerate(key):
		decimal += int(bit)*2**(-1-i)
	if value != 0:
		res_decimal[str(decimal)] += value
		res_energy[str(-decimal*2*np.pi/dt)] += value

print(res_decimal)
print(res_energy)


