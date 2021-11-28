import velociraptor as vr
import unyt
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
data = vr.load('/data3/FLAMINGO/REF/L0200N0360/VR/halos_0008.properties.0')

halo_masses = data.masses.mass_200crit
halo_xpositions = data.positions.xc
halo_ypositions = data.positions.yc
halo_zpositions = data.positions.zc

a = 0
N1 = 0
N2 = 2000
N = N1+N2
Boxsize1 = 200
Boxsize2 = 20

def distance(v1,v2,boxsize): 
	distance =(((v1[0] - v2[0] + boxsize / 2) % boxsize - boxsize / 2) ** 2 + ((v1[1] - v2[1] + boxsize / 2) % boxsize - boxsize / 2) ** 2 + ((v1[2] - v2[2] + boxsize / 2) % boxsize - boxsize / 2) ** 2) ** 0.5
	return distance

def nonpdistance(v1,v2):
	distance =((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2) ** 0.5
	return distance

position = []

for i in range (0, len(data.ids.id) - 1):
	if data.structure_type.structuretype[i] == 10:
		if 0.9e+12 < halo_masses.to_value(unyt.msun)[i] < 1.1e+12:
			a = a+1
			position.append([halo_xpositions.to_value(unyt.Mpc)[i], halo_ypositions.to_value(unyt.Mpc)[i], halo_zpositions.to_value(unyt.Mpc)[i]])
			if a == 1000:
				break

print(a)

D = position
R = np.random.uniform(0,20,(N2,3))

rdd = []
rrr = []
rdr = []
b = 0
c = 0

for i in range(a):
	for j in range(i+1,a):
			r1 = distance(D[i],D[j],Boxsize1)
			if r1 <= 20:
				b += 1
				rdd.append(r1)
print(b)
for i in range(N):
	for j in range(i+1,N):
			r2 = distance(R[i], R[j], Boxsize2)
			if r2 <= 20:
				c += 1
				rrr.append(r2)
print(c)
'''
for i in range(N):
	for j in range(a):
		r3 = distance(R[i], D[j])
		rdr.append(r3)
'''
histdd, bin_edge = np.histogram(rdd, bins=np.arange(1,20,0.5))
histrr, bin_edge = np.histogram(rrr, bins=np.arange(1,20,0.5))
# histdr, bin_edge = np.histogram(rdr, bins=np.arange(1,30,0.5))

# estimator = [(((histdd[i]/(a*(a-1))) - (histdr[i]/(N*a)))/(histrr[i]/(N*(N-1)))) + 1 for i in range(0,len(histdd))]
estimator = [((histdd[i]/(a*(a-1)/2))/(histrr[i]/(2000*(2000000 -1)/2))) - 1 for i in range(0,len(histdd))]

X = (bin_edge[:-1] + bin_edge[1:])/2
Y = estimator
Y1 = (X/5) ** (-1.7)
print(X)
print(Y)

plt.scatter(X,Y, label='0.9*10^12<M_halo<1.1*10^12', c='b')
plt.plot(X,Y)
plt.plot(X,Y1, label='good fit $\zeta$=r^(-1.7)')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-2,1e2)
#plt.xlim(1e-1,1e2)
plt.xlabel('log(r)')
plt.ylabel(r'log($\zeta$)')
plt.legend()
plt.savefig('correlation_halo_log_N_4000_2000_1000_periodic_arange_uniform_simpleestimator_changeboxsize_smallboxperiodic20.png')
