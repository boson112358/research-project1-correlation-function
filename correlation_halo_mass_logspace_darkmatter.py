import velociraptor as vr
import unyt
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
data = vr.load('/data3/FLAMINGO/DMO/L0200N0360/VR/halos_0008.properties.0')

halo_masses = data.masses.mass_200crit
halo_xpositions = data.positions.xc
halo_ypositions = data.positions.yc
halo_zpositions = data.positions.zc

a = 0
b = 0
c = 0
N = 3000
Boxsize = 200 

def distance(v1,v2): 
	distance =(((v1[0] - v2[0] + Boxsize / 2) % Boxsize - Boxsize / 2) ** 2 + ((v1[1] - v2[1] + Boxsize / 2) % Boxsize - Boxsize / 2) ** 2 + ((v1[2] - v2[2] + Boxsize / 2) % Boxsize - Boxsize / 2) ** 2) ** 0.5
	return distance

position1 = []
position2 = []
position3 = []

for i in range (0, len(data.ids.id)):
	if data.structure_type.structuretype[i] == 10:
		if 0.9e+11 < halo_masses.to_value(unyt.msun)[i] < 1.1e+11:
			a = a+1
			position1.append([halo_xpositions.to_value(unyt.Mpc)[i], halo_ypositions.to_value(unyt.Mpc)[i], halo_zpositions.to_value(unyt.Mpc)[i]])
			if a == 1000:
				break
print(a)

for i in range(0, len(data.ids.id)):
	if data.structure_type.structuretype[i] == 10:
		if 0.9e+12 < halo_masses.to_value(unyt.msun)[i] < 1.1e+12:
			b = b+1
			position2.append([halo_xpositions.to_value(unyt.Mpc)[i], halo_ypositions.to_value(unyt.Mpc)[i], halo_zpositions.to_value(unyt.Mpc)[i]])
			if b == 1000:
				break
print(b)

for i in range(0, len(data.ids.id)):
	if data.structure_type.structuretype[i] == 10:
		if 0.6e+13 < halo_masses.to_value(unyt.msun)[i] < 1.4e+13:
			c = c+1
			position3.append([halo_xpositions.to_value(unyt.Mpc)[i], halo_ypositions.to_value(unyt.Mpc)[i], halo_zpositions.to_value(unyt.Mpc)[i]])
			if c == 1000:
				break
print(c)

D1 = position1
D2 = position2
D3 = position3
R = np.random.uniform(0,200,(N,3))

matrixdd1 = np.zeros((a,a))
matrixdd2 = np.zeros((b,b))
matrixdd3 = np.zeros((c,c))
matrixrr = np.zeros((N,N))
matrixdr1 = np.zeros((a,N))
matrixdr2 = np.zeros((b,N))
matrixdr3 = np.zeros((c,N))


for i in range(a):
	for j in range(i+1,a):
		matrixdd1[i][j] = distance(D1[i], D1[j])

for i in range(b):
	for j in range(i+1,b):
		matrixdd2[i][j] = distance(D2[i], D2[j])

for i in range(c):
	for j in range(i+1,c):
		matrixdd3[i][j] = distance(D3[i], D3[j])

for i in range(N):
	for j in range(i+1,N):
		matrixrr[i][j] = distance(R[i], R[j])

for i in range(a):
	for j in range(N):
		matrixdr1[i][j] = distance(D1[i], R[j])

for i in range(b):
	for j in range(N):
		matrixdr2[i][j] = distance(D2[i], R[j])

for i in range(c):
	for j in range(N):
		matrixdr3[i][j] = distance(D3[i], R[j])



histdd1, bin_edge = np.histogram(matrixdd1, bins=np.arange(1,30,0.5))
histdd2, bin_edge = np.histogram(matrixdd2, bins=np.arange(1,30,0.5))
histdd3, bin_edge = np.histogram(matrixdd3, bins=np.arange(1,30,0.5))
histrr, bin_edge = np.histogram(matrixrr, bins=np.arange(1,30,0.5))
histdr1, bin_edge = np.histogram(matrixdr1, bins=np.arange(1,30,0.5))
histdr2, bin_edge = np.histogram(matrixdr2, bins=np.arange(1,30,0.5))
histdr3, bin_edge = np.histogram(matrixdr3, bins=np.arange(1,30,0.5))

estimator1 = [(((histdd1[i]/(a*(a-1))) - (histdr1[i]/(N*a)))/(histrr[i]/(N*(N-1)))) + 1 for i in range(0,len(histdd1))]
estimator2 = [(((histdd2[i]/(b*(b-1))) - (histdr2[i]/(N*b)))/(histrr[i]/(N*(N-1)))) + 1 for i in range(0,len(histdd1))]
estimator3 = [(((histdd3[i]/(c*(c-1))) - (histdr3[i]/(N*c)))/(histrr[i]/(N*(N-1)))) + 1 for i in range(0,len(histdd1))]


X = (bin_edge[:-1] + bin_edge[1:])/2
Y1 = estimator1
Y2 = estimator2
Y3 = estimator3
print(Y1)
print(Y2)
print(Y3)

plt.scatter(X,Y1, label='0.9*10^11<M_halo<1.1*10^11')
plt.scatter(X, Y2, label='0.9*10^12<M_halo<1.1*10^12')
plt.scatter(X, Y3, label='0.6*10^13<M_halo<1.4*10^13')
plt.plot(X,Y1,X,Y2,X,Y3)
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-2, 1e2)
plt.xlabel('log(r)')
plt.ylabel(r'log($\zeta$)')
plt.legend()
plt.savefig('correlation_halo_Mass3_N_3000_1000_periodic_logspace_darkmatter_1.png')
