import qiskit as qk
import matplotlib.pylab as plt
import numpy as np
qk.IBMQ.load_accounts()

class QCPairing:
	def __init__(self,n_work,n_simulation,delta=1,g=1,dt=0.001,Emax=500):
		self.n_work = n_work
		self.n_simulation = n_simulation
		self.n_qubits = n_work + n_simulation + 1
		self.qb = qk.QuantumRegister(self.n_qubits)
		self.cb = qk.ClassicalRegister(self.n_qubits)
		self.qz = qk.QuantumCircuit(self.qb,self.cb)
		self.delta = delta
		self.g = g
		self.dt = dt
		self.Emax = Emax
	def set_dt(self,dt):
		self.dt = dt
	def set_Emax(self,Emax):
		self.Emax = Emax
	def set_delta(self,delta):
		self.delta = delta
	def set_g(self,g):
		self.g = g
	
	def H0(self,dt,control_qubit):
		g = self.g
		delta=self.delta
		n_work = self.n_work
		n_simulation=self.n_simulation
		n_qubits=self.n_qubits
		Emax=self.Emax
		qb = self.qb
		cb = self.cb
		qz = self.qz
		s_state = 0
		for q_state in range(0,n_simulation):
			if q_state % 2 == 0:
				s_state += 1
			qz.crz(dt*delta*(s_state - 1),qb[control_qubit],qb[q_state+n_work])

		qz.cu1(-dt*delta*(1/8)*(n_simulation-2)*n_simulation,qb[control_qubit],qb[n_work])
		qz.x(qb[n_work])
		qz.cu1(-dt*delta*(1/8)*(n_simulation-2)*n_simulation,qb[control_qubit],qb[n_work])
		qz.x(qb[n_work])

		qz.cu1(Emax*dt,qb[control_qubit],qb[n_work])
		qz.x(qb[n_work])
		qz.cu1(Emax*dt,qb[control_qubit],qb[n_work])
		qz.x(qb[n_work])

		self.qb = qb
		self.cb = cb
		self.qz = qz

	def H1(self,dt,control_qubit):
		g = self.g
		delta=self.delta
		n_work = self.n_work
		n_simulation=self.n_simulation
		n_qubits=self.n_qubits
		qb = self.qb
		cb = self.cb
		qz = self.qz
		for p in range(1,n_simulation,2):
			for q in range(p,n_simulation,2):
				if p == q:
					theta = -2*(1/8)*g*dt
					qz.cu1(-theta/2,qb[control_qubit],qb[n_work])
					qz.x(qb[n_work])
					qz.cu1(-theta/2,qb[control_qubit],qb[n_work])
					qz.x(qb[n_work])

					qz.crz(theta,qb[control_qubit],qb[p-1+n_work])
					qz.crz(theta,qb[control_qubit],qb[p+n_work])
					
					qz.cx(qb[p+n_work],qb[n_qubits-1])
					qz.cx(qb[p-1+n_work],qb[n_qubits-1])
					qz.crz(-theta,qb[control_qubit],qb[n_qubits-1])
					qz.cx(qb[p+n_work],qb[n_qubits-1])
					qz.cx(qb[p-1+n_work],qb[n_qubits-1])
				else:
					theta = -2*(1/16)*g*dt
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
		self.qb = qb
		self.cb = cb
		self.qz = qz


	def PhaseEstimation(self,t):
		dt = self.dt
		g = self.g
		delta=self.delta
		n_work = self.n_work
		n_simulation=self.n_simulation
		n_qubits=self.n_qubits
		qb = self.qb
		cb = self.cb
		qz = self.qz

		#Initialize simulation qubits to superposition of bell states

		for i in range(0,n_simulation,2):
			qz.h(qb[n_work+i])
			qz.cx(qb[n_work+i],qb[n_work+i+1])
		



		for cq in range(n_work):
			qz.h(qb[cq])


		for cq in range(n_work):
			for j in range(int(t/dt)):
				self.H0((2**cq)*dt,cq)
				self.H1((2**cq)*dt,cq)

		self.qb = qb
		self.cb = cb
		self.qz = qz

	def inverse_Fourier(self):
		dt = self.dt
		g = self.g
		delta=self.delta
		n_work = self.n_work
		n_simulation=self.n_simulation
		n_qubits=self.n_qubits
		qb = self.qb
		cb = self.cb
		qz = self.qz
		for cq in range(int(n_work/2)):
			qz.swap(qb[cq],qb[n_work-cq-1])
		for cq in range(n_work):
			for i in range(cq):
				qz.cu1(-2*np.pi/(2**(1+cq-i)),qb[i],qb[cq])
			qz.h(qb[cq])
		self.qb = qb
		self.cb = cb
		self.qz = qz

	def solve(self,t=0):
		if t == 0:
			t = self.dt
		self.PhaseEstimation(t)
		self.inverse_Fourier()
		return(self.qz,self.qb,self.cb)


#[-1.11803399  1.11803399]
n_work = 8
n_simulation = 4
Emax=5
dt = 0.001
g = 1
delta=1
t= 100*dt
Pairing1 = QCPairing(n_work,n_simulation,delta=delta,g=g,dt=dt,Emax=Emax)
qz,qb,cb= Pairing1.solve(t)




qz.measure(qb,cb)


shots = 1000
job = qk.execute(qz, backend = qk.Aer.get_backend('qasm_simulator'), shots=shots)
result = job.result()

result = result.get_counts(qz)
res_decimal = {}
res_energy = {}
for key,value in result.items():
	key = key[n_simulation+1:]
	decimal = 0
	for i,bit in enumerate(key):
		decimal += int(bit)*2**(-i-1)
	if value != 0:
		res_decimal[decimal] = 0
		res_energy[Emax-decimal*2*np.pi/t] = 0
for key,value in result.items():
	key = key[n_simulation+1:]
	decimal = 0
	for i,bit in enumerate(key):
		decimal += int(bit)*2**(-1-i)
	if value != 0:
		res_decimal[decimal] += value
		res_energy[Emax-decimal*2*np.pi/t] += value

print(res_decimal)
print(res_energy)

lists = sorted(res_energy.items()) # sorted by key, return a list of tuples

x, y = zip(*lists)

plt.plot(x,y)
plt.xlabel('Eigenvalue')
plt.ylabel('Times measured')
plt.show()

