import velociraptor as vr
import unyt
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import log
data1 = vr.load('/data3/FLAMINGO/REF/L0200N0360/VR/halos_0008.properties.0')
data2 = vr.load('/data3/FLAMINGO/DMO/L0200N0360/VR/halos_0008.properties.0')

halo_masses1 = data1.masses.mass_200crit
halo_xpositions1 = data1.positions.xc
halo_ypositions1 = data1.positions.yc
halo_zpositions1 = data1.positions.zc

halo_masses2 = data2.masses.mass_200crit
halo_xpositions2 = data2.positions.xc
halo_ypositions2 = data2.positions.yc
halo_zpositions2 = data2.positions.zc

a = 0
b = 0
N1 = 0
N2 = 1000
N = N1+N2
Boxsize1 = 200
Boxsize2 = 20

def distance(v1,v2,boxsize): 
	distance =(((v1[0] - v2[0] + boxsize / 2) % boxsize - boxsize / 2) ** 2 + ((v1[1] - v2[1] + boxsize / 2) % boxsize - boxsize / 2) ** 2 + ((v1[2] - v2[2] + boxsize / 2) % boxsize - boxsize / 2) ** 2) ** 0.5
	return distance

def nonpdistance(v1,v2):
	distance =((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2) ** 0.5
	return distance
'''
def fund(x,a,b):
	return x**a + b
'''
def fund(x,a,b):
	return a*log(x) + b


position1 = []
position2 = []

# import halo position
for i in range (0, len(data1.ids.id) - 1):
	if data1.structure_type.structuretype[i] == 10:
		if 0.9e+12 < halo_masses1.to_value(unyt.msun)[i] < 1.1e+12:
			a = a+1
			position1.append([halo_xpositions1.to_value(unyt.Mpc)[i], halo_ypositions1.to_value(unyt.Mpc)[i], halo_zpositions1.to_value(unyt.Mpc)[i]])
			if a == 2000:
				break
print(a)

for i in range (0, len(data2.ids.id) - 1):
	if data2.structure_type.structuretype[i] == 10:
		if 0.9e+12 < halo_masses2.to_value(unyt.msun)[i] < 1.1e+12:
			b = b+1
			position2.append([halo_xpositions2.to_value(unyt.Mpc)[i], halo_ypositions2.to_value(unyt.Mpc)[i], halo_zpositions2.to_value(unyt.Mpc)[i]])
			if b == 2000:
				break



D1 = position1
D2 = position2
R = np.random.uniform(0,20,(N2,3))

matrixdd1 = np.zeros((a,a))
matrixdd2 = np.zeros((b,b))
matrixrr = np.zeros((N2,N2))
# matrixdr = np.zeros((a,N2))

for i in range(a):
	for j in range(i+1,a):
		r1 = distance(D1[i],D1[j],Boxsize1)
		if r1 <= 20:
			matrixdd1[i][j] = r1
for i in range(b):
	for j in range(i+1,b):
		r2 = distance(D2[i],D2[j],Boxsize1)
		if r2 <= 20:
			matrixdd2[i][j] = r2


for i in range(N):
	for j in range(i+1,N):
			r2 = distance(R[i], R[j], Boxsize2)
			if r2 <= 20:
				matrixrr[i][j] = r2

'''
for i in range(N):
	for j in range(a):
		r3 = distance(R[i], D[j])
		rdr.append(r3)
'''
histdd1, bin_edge = np.histogram(matrixdd1, bins=np.arange(1,10,0.2))
histdd2, bin_edge = np.histogram(matrixdd2, bins=np.arange(1,10,0.2))
histrr, bin_edge = np.histogram(matrixrr, bins=np.arange(1,10,0.2))
# histdr, bin_edge = np.histogram(rdr, bins=np.arange(1,30,0.5))

# estimator = [(((histdd[i]/(a*(a-1))) - (histdr[i]/(N*a)))/(histrr[i]/(N*(N-1)))) + 1 for i in range(0,len(histdd))]
# simple estimator
estimator1 = [((histdd1[i]/(a*(a-1)/2))/(histrr[i]/(N2*(N2*1000 -1)/2))) - 1 for i in range(0,len(histdd1))]
estimator2 = [((histdd2[i]/(a*(a-1)/2))/(histrr[i]/(N2*(N2*1000 -1)/2))) - 1 for i in range(0,len(histdd2))]

X = (bin_edge[:-1] + bin_edge[1:])/2
Y1 = estimator1
Y2 = estimator2
print(X)
print(Y1)
print(Y2)

# fit curve
popt1, pcov = curve_fit(fund, X,log(Y1))
popt2, pcov = curve_fit(fund, X,log(Y2))
Yf1 = [math.exp(fund(i, popt1[0], popt1[1])) for i in X]
Yf2 = [math.exp(fund(i, popt2[0], popt2[1])) for i in X]
# Y1 = (X/5) ** (-1.7)

plt.scatter(X,Y1, label='0.9*10^12<M_halo<1.1*10^12', c='b')
plt.scatter(X,Y2, label='DMO 0.9*10^12<M_halo<1.1*10^12', c='y')
plt.plot(X,Y1,X,Y2)
plt.plot(X,Yf1, label='fit curve log($\zeta$) = %5.3f*log(r) + %5.3f'% tuple(popt1))
plt.plot(X,Yf2, label='fit curve log($\zeta$) = %5.3f*log(r) + %5.3f'% tuple(popt2))
# plt.plot(X,Y1, label='good fit $\zeta$=r^(-1.7)')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-2,1e2)
#plt.xlim(1e-1,1e2)
plt.xlabel('log(r)')
plt.ylabel(r'log($\zeta$)')
plt.legend()
plt.savefig('correlation_halo_log_N_2000_1000_periodic_arange_uniform_simpleestimator_smallboxperiodic20_fitcurve_compareDMO.png')
