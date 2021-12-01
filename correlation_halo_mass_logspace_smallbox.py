import velociraptor as vr
import unyt
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import log
data = vr.load('/data3/FLAMINGO/REF/L0200N0360/VR/halos_0008.properties.0')

halo_masses = data.masses.mass_200crit
halo_xpositions = data.positions.xc
halo_ypositions = data.positions.yc
halo_zpositions = data.positions.zc

a = 0
b = 0
c = 0
d = 0
e = 0
f = 0
N1 = 0
N2 = 1000
N = N1+N2
Boxsize1 = 200
Boxsize2 = 20

def distance(v1,v2,boxsize): 
	distance =(((v1[0] - v2[0] + boxsize / 2) % boxsize - boxsize / 2) ** 2 + ((v1[1] - v2[1] + boxsize / 2) % boxsize - boxsize / 2) ** 2 + ((v1[2] - v2[2] + boxsize / 2) % boxsize - boxsize / 2) ** 2) ** 0.5
	return distance

def fund(x,a,b):
	return a*log(x) + b

position1 = []
position2 = []
position3 = []
position4 = []
position5 = []
position6 = []

for i in range (0, len(data.ids.id)):
	if data.structure_type.structuretype[i] == 10:
		if 1e+11 < halo_masses.to_value(unyt.msun)[i] < 3.16e+11:
			a = a+1
			position1.append([halo_xpositions.to_value(unyt.Mpc)[i], halo_ypositions.to_value(unyt.Mpc)[i], halo_zpositions.to_value(unyt.Mpc)[i]])
			if a == 2000:
				break
print(a)

for i in range(0, len(data.ids.id)):
	if data.structure_type.structuretype[i] == 10:
		if 3.16e+11 < halo_masses.to_value(unyt.msun)[i] < 1e+12:
			b = b+1
			position2.append([halo_xpositions.to_value(unyt.Mpc)[i], halo_ypositions.to_value(unyt.Mpc)[i], halo_zpositions.to_value(unyt.Mpc)[i]])
			if b == 2000:
				break
print(b)

for i in range(0, len(data.ids.id)):
	if data.structure_type.structuretype[i] == 10:
		if 1e+12 < halo_masses.to_value(unyt.msun)[i] < 3.16e+12:
			c = c+1
			position3.append([halo_xpositions.to_value(unyt.Mpc)[i], halo_ypositions.to_value(unyt.Mpc)[i], halo_zpositions.to_value(unyt.Mpc)[i]])
			if c == 2000:
				break
print(c)

for i in range(0, len(data.ids.id)):
	if data.structure_type.structuretype[i] == 10:
		if 3.16e+12 < halo_masses.to_value(unyt.msun)[i] < 1e+13:
			d = d+1
			position4.append([halo_xpositions.to_value(unyt.Mpc)[i], halo_ypositions.to_value(unyt.Mpc)[i], halo_zpositions.to_value(unyt.Mpc)[i]])
			if d == 2000:
				break
print(d)
'''
for i in range(0, len(data.ids.id)):
	if data.structure_type.structuretype[i] == 10:
		if 1e+13 < halo_masses.to_value(unyt.msun)[i] < 3.16e+13:
			e = e+1
			position5.append([halo_xpositions.to_value(unyt.Mpc)[i], halo_ypositions.to_value(unyt.Mpc)[i], halo_zpositions.to_value(unyt.Mpc)[i]])
			if e == 1000:
				break
print(e)

for i in range(0, len(data.ids.id)):
	if data.structure_type.structuretype[i] == 10:
		if 3.16e+13 < halo_masses.to_value(unyt.msun)[i] < 1e+14:
			f = f+1
			position6.append([halo_xpositions.to_value(unyt.Mpc)[i], halo_ypositions.to_value(unyt.Mpc)[i], halo_zpositions.to_value(unyt.Mpc)[i]])
			if f == 1000:
				break
print(f)
'''
D1 = position1
D2 = position2
D3 = position3
D4 = position4
D5 = position5
D6 = position6
R = np.random.uniform(0,20,(N2,3))

matrixdd1 = np.zeros((a,a))
matrixdd2 = np.zeros((b,b))
matrixdd3 = np.zeros((c,c))
matrixdd4 = np.zeros((d,d))
matrixdd5 = np.zeros((e,e))
matrixdd6 = np.zeros((f,f))
matrixrr = np.zeros((N,N))
matrixdr1 = np.zeros((a,N))
matrixdr2 = np.zeros((b,N))
matrixdr3 = np.zeros((c,N))


for i in range(a):
	for j in range(i+1,a):
		if distance(D1[i], D1[j], Boxsize1) <= 20:
			matrixdd1[i][j] = distance(D1[i], D1[j], Boxsize1)

for i in range(b):
	for j in range(i+1,b):
		if distance(D2[i], D2[j], Boxsize1) <= 20:
			matrixdd2[i][j] = distance(D2[i], D2[j], Boxsize1)

for i in range(c):
	for j in range(i+1,c):
		if distance(D3[i],D3[j], Boxsize1) <= 20:
			matrixdd3[i][j] = distance(D3[i], D3[j], Boxsize1)

for i in range(d):
	for j in range(i+1,d):
		if distance(D4[i],D4[j], Boxsize1) <= 20:
			matrixdd4[i][j] = distance(D4[i], D4[j], Boxsize1)
'''
for i in range(e):
	for j in range(i+1,e):
		if distance(D5[i],D5[j], Boxsize1) <= 20:
			matrixdd5[i][j] = distance(D5[i], D5[j], Boxsize1)

for i in range(f):
	for j in range(i+1,f):
		if distance(D6[i],D6[j], Boxsize1) <= 20:
			matrixdd6[i][j] = distance(D6[i], D6[j], Boxsize1)
'''

for i in range(N):
	for j in range(i+1,N):
		if distance(R[i], R[j], Boxsize2):
			matrixrr[i][j] = distance(R[i], R[j], Boxsize2)
'''
for i in range(a):
	for j in range(N):
		matrixdr1[i][j] = distance(D1[i], R[j])

for i in range(b):
	for j in range(N):
		matrixdr2[i][j] = distance(D2[i], R[j])

for i in range(c):
	for j in range(N):
		matrixdr3[i][j] = distance(D3[i], R[j])
'''


histdd1, bin_edge = np.histogram(matrixdd1, bins=np.arange(1,10,0.2))
histdd2, bin_edge = np.histogram(matrixdd2, bins=np.arange(1,10,0.2))
histdd3, bin_edge = np.histogram(matrixdd3, bins=np.arange(1,10,0.2))
histdd4, bin_edge = np.histogram(matrixdd4, bins=np.arange(1,10,0.2))
#histdd5, bin_edge = np.histogram(matrixdd5, bins=np.arange(1,10,0.2))
#histdd6, bin_edge = np.histogram(matrixdd6, bins=np.arange(1,10,0.2))
histrr, bin_edge = np.histogram(matrixrr, bins=np.arange(1,10,0.2))
#histdr1, bin_edge = np.histogram(matrixdr1, bins=np.arange(1,30,0.5))
#histdr2, bin_edge = np.histogram(matrixdr2, bins=np.arange(1,30,0.5))
#histdr3, bin_edge = np.histogram(matrixdr3, bins=np.arange(1,30,0.5))

#estimator1 = [(((histdd1[i]/(a*(a-1))) - (histdr1[i]/(N*a)))/(histrr[i]/(N*(N-1)))) + 1 for i in range(0,len(histdd1))]
#estimator2 = [(((histdd2[i]/(b*(b-1))) - (histdr2[i]/(N*b)))/(histrr[i]/(N*(N-1)))) + 1 for i in range(0,len(histdd1))]
#estimator3 = [(((histdd3[i]/(c*(c-1))) - (histdr3[i]/(N*c)))/(histrr[i]/(N*(N-1)))) + 1 for i in range(0,len(histdd1))]

estimator1 = [((histdd1[i]/(a*(a-1)/2))/(histrr[i]/(1000*(1000000-1)/2))) - 1 for i in range(0,len(histdd1))]
estimator2 = [((histdd2[i]/(a*(a-1)/2))/(histrr[i]/(1000*(1000000-1)/2))) - 1 for i in range(0,len(histdd1))]
estimator3 = [((histdd3[i]/(a*(a-1)/2))/(histrr[i]/(1000*(1000000-1)/2))) - 1 for i in range(0,len(histdd1))]
estimator4 = [((histdd4[i]/(a*(a-1)/2))/(histrr[i]/(1000*(1000000-1)/2))) - 1 for i in range(0,len(histdd1))]
#estimator5 = [((histdd3[i]/(a*(a-1)/2))/(histrr[i]/(1000*(1000000-1)/2))) - 1 for i in range(0,len(histdd1))]
#estimator6 = [((histdd3[i]/(a*(a-1)/2))/(histrr[i]/(1000*(1000000-1)/2))) - 1 for i in range(0,len(histdd1))]



X = (bin_edge[:-1] + bin_edge[1:])/2
Y1 = estimator1
Y2 = estimator2
Y3 = estimator3
Y4 = estimator4
#Y5 = estimator5
#Y6 = estimator6
print(Y1)
print(Y2)
print(Y3)
print(Y4)
#print(Y5)
#print(Y6)

# Y4 = (X/5) ** (-1.7)
# fit curve
popt1, pcov1 = curve_fit(fund,X,log(Y1))
popt2, pcov1 = curve_fit(fund,X,log(Y2))
popt3, pcov1 = curve_fit(fund,X,log(Y3))
popt4, pcov1 = curve_fit(fund,X,log(Y4))
#popt5, pcov1 = curve_fit(fund,X,log(Y5))
#popt6, pcov1 = curve_fit(fund,X,log(Y6))

Y1f = [math.exp(fund(i, popt1[0], popt1[1])) for i in X]
Y2f = [math.exp(fund(i, popt2[0], popt2[1])) for i in X]
Y3f = [math.exp(fund(i, popt3[0], popt3[1])) for i in X]
Y4f = [math.exp(fund(i, popt4[0], popt4[1])) for i in X]
#Y5f = [math.exp(fund(i, popt5[0], popt5[1])) for i in X]
#Y6f = [math.exp(fund(i, popt6[0], popt6[1])) for i in X]

plt.scatter(X,Y1, label='11.0 < log(${M_{halo}}$) < 11.5')
plt.scatter(X, Y2, label='11.5 < log(${M_{halo}}$) < 12')
plt.scatter(X, Y3, label='12.0 < log(${M_{halo}}$) < 12.5')
plt.scatter(X, Y4, label='12.5 < log(${M_{halo}}$) < 13.0')
#plt.scatter(X, Y5, label='13.0 < log(${M_{halo}}$) < 13.5')
#plt.scatter(X, Y6, label='13.5 < log(${M_{halo}}$) < 14.0')
plt.plot(X,Y1,X,Y2,X,Y3,X,Y4)
plt.plot(X,Y1f, label='fit curve1 log($\zeta$)=%5.3flog(r) + %5.3f'%tuple(popt1))
plt.plot(X,Y2f, label='fit curve2 log($\zeta$)=%5.3flog(r) + %5.3f'%tuple(popt2))
plt.plot(X,Y3f, label='fit curve3 log($\zeta$)=%5.3flog(r) + %5.3f'%tuple(popt3))
plt.plot(X,Y4f, label='fit curve4 log($\zeta$)=%5.3flog(r) + %5.3f'%tuple(popt4))
#plt.plot(X,Y5f, label='fit curve5 log($\zeta$)=%5.3flog(r) + %5.3f'%tuple(popt5))
#plt.plot(X,Y6f, label='fit curve6 log($\zeta$)=%5.3flog(r) + %5.3f'%tuple(popt6))
# plt.plot(X,Y4, label='good fit $\zeta$=(r/5)^(-1.7)')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-2,1e2)
plt.xlabel('log(r)')
plt.ylabel(r'log($\zeta$)')
plt.legend()
plt.savefig('correlation_halo_Mass3_N_2000_1000_periodic_arrange_smallbox_fitcurve_6.png')
