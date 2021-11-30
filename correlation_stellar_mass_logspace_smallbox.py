import velociraptor as vr
import unyt
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
data = vr.load('/data3/FLAMINGO/REF/L0200N0360/VR/halos_0008.properties.0')

stellar_masses = data.apertures.mass_star_30_kpc
halo_xpositions = data.positions.xc
halo_ypositions = data.positions.yc
halo_zpositions = data.positions.zc

a = 0
b = 0
c = 0
N2 = 1000
Boxsize1 = 200 
Boxsize2 = 20

def distance(v1,v2,boxsize): 
	distance =(((v1[0] - v2[0] + boxsize / 2) % boxsize - boxsize / 2) ** 2 + ((v1[1] - v2[1] + boxsize / 2) % boxsize - boxsize / 2) ** 2 + ((v1[2] - v2[2] + boxsize / 2) % boxsize - boxsize / 2) ** 2) ** 0.5
	return distance

position1 = []
position2 = []
position3 = []


for i in range(0, len(data.ids.id)):
	if 1e+10 < stellar_masses.to_value(unyt.msun)[i] < 3.16e+10:
		a = a+1
		position1.append([halo_xpositions.to_value(unyt.Mpc)[i], halo_ypositions.to_value(unyt.Mpc)[i], halo_zpositions.to_value(unyt.Mpc)[i]])
		if a == 2000:
			break
print(a)

for i in range(0, len(data.ids.id)):
	if 3.16e+10 < stellar_masses.to_value(unyt.msun)[i] < 1e+11:
		b = b+1
		position2.append([halo_xpositions.to_value(unyt.Mpc)[i], halo_ypositions.to_value(unyt.Mpc)[i], halo_zpositions.to_value(unyt.Mpc)[i]])
		if b == 2000:
			break
print(b)

for i in range(0, len(data.ids.id)):
	if 1e+11 < stellar_masses.to_value(unyt.msun)[i] < 3.16e+11:
		c = c+1
		position3.append([halo_xpositions.to_value(unyt.Mpc)[i], halo_ypositions.to_value(unyt.Mpc)[i], halo_zpositions.to_value(unyt.Mpc)[i]])
		if c == 2000:
			break
print(c)



D1 = position1
D2 = position2
D3 = position3
R = np.random.uniform(0,20,(N2,3))

matrixdd1 = np.zeros((a,a))
matrixdd2 = np.zeros((b,b))
matrixdd3 = np.zeros((c,c))
matrixrr = np.zeros((N2,N2))
#matrixdr1 = np.zeros((a,N))
#matrixdr2 = np.zeros((b,N))
#matrixdr3 = np.zeros((c,N))


for i in range(a):
	for j in range(i+1,a):
		if distance(D1[i],D1[j],Boxsize1) <= 20:
			matrixdd1[i][j] = distance(D1[i], D1[j],Boxsize1)

for i in range(b):
	for j in range(i+1,b):
		if distance(D2[i],D2[j],Boxsize1) <= 20:
			matrixdd2[i][j] = distance(D2[i], D2[j],Boxsize1)

for i in range(c):
	for j in range(i+1,c):
		if distance(D3[i], D3[j],Boxsize1) <= 20:
			matrixdd3[i][j] = distance(D3[i], D3[j],Boxsize1)

for i in range(N2):
	for j in range(i+1,N2):
		if distance(R[i],R[j],Boxsize2) <= 20:
			matrixrr[i][j] = distance(R[i], R[j],Boxsize2)
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


histdd1, bin_edge = np.histogram(matrixdd1, bins=np.arange(0.1,10,0.2))
histdd2, bin_edge = np.histogram(matrixdd2, bins=np.arange(0.1,10,0.2))
histdd3, bin_edge = np.histogram(matrixdd3, bins=np.arange(0.1,10,0.2))
histrr, bin_edge = np.histogram(matrixrr, bins=np.arange(0.1,10,0.2))
#histdr1, bin_edge = np.histogram(matrixdr1, bins=np.logspace(0,1.5,num=30))
#histdr2, bin_edge = np.histogram(matrixdr2, bins=np.arange(1,30,0.5))
#histdr3, bin_edge = np.histogram(matrixdr3, bins=np.arange(1,30,0.5))

#estimator1 = [(((histdd1[i]/(a*(a-1))) - (histdr1[i]/(N*a)))/(histrr[i]/(N*(N-1)))) + 1 for i in range(0,len(histdd1))]
#estimator2 = [(((histdd2[i]/(b*(b-1))) - (histdr2[i]/(N*b)))/(histrr[i]/(N*(N-1)))) + 1 for i in range(0,len(histdd2))]
#estimator3 = [(((histdd3[i]/(c*(c-1))) - (histdr3[i]/(N*c)))/(histrr[i]/(N*(N-1)))) + 1 for i in range(0,len(histdd2))]

estimator1 = [((histdd1[i]/(a*(a-1)/2))/(histrr[i]/(1000*(1000000-1)/2))) - 1 for i in range(0,len(histdd2))]
estimator2 = [((histdd2[i]/(b*(b-1)/2))/(histrr[i]/(1000*(1000000-1)/2))) - 1 for i in range(0,len(histdd2))]
estimator3 = [((histdd3[i]/(c*(c-1)/2))/(histrr[i]/(1000*(1000000-1)/2))) - 1 for i in range(0,len(histdd2))]

X = (bin_edge[:-1] + bin_edge[1:])/2
Y1 = estimator1
Y2 = estimator2
Y3 = estimator3
print(Y1)
print(Y2)
print(Y3)

# Y4 = (X/5) ** (-1.7)

plt.scatter(X,Y1, label='10.0 < log(${M_*}$(30kpc)) < 10.5')
plt.scatter(X, Y2, label='10.5 < log(${M_*}$(30kpc)) < 11')
plt.scatter(X, Y3, label='11 < log(${M_*}$(30kpc)) < 11.5')
plt.plot(X,Y1,X,Y2,X,Y3)
# plt.plot(X,Y4,label='good fit $\zeta$=(r/5)^(-1.7)')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-2, 1e3)
plt.xlabel('log(r)')
plt.ylabel(r'log($\zeta$)')
plt.legend()
plt.savefig('correlation_stellar_Mass3_N_2000_1000_periodic_smallbox_below1_3.png')
