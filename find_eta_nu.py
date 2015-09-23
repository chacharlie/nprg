from pylab import *
from numpy import *
import  matplotlib.pyplot as plt

model=input('Enter model (1 for ON, 2 for A):')
if model==1:
	from load_results import *
elif model==2:
	from A_load_results import *

maxEtaZ=[]
maxEtaX=[]
maxNu=[]

for i in range(Ndicho):
	EtaZt=matrixEtaZ[i]
	EtaXt=matrixEtaX[i]
	#dtEtaZ = abs(gradient(EtaZt)/dt)
	#ind = argmin(dtEtaZ)
	#maxEtaZ.append(EtaZt[ind])
	maxEtaZ.append(EtaZt[argmin(abs(gradient(EtaZt)))])
	maxEtaX.append(EtaXt[argmin(abs(gradient(EtaXt)))])


indexVj=0
for istep in range(Ndicho):
	Vj=[]
  	NTi=len(matrixEtaZ[istep])  
  	ti=linspace(0,100*dt*NTi,NTi)
	for k in range(len(matrixV[istep])):
		Vj.append(matrixV[istep][k][indexVj].real)
	dlndVj=gradient(log(abs(gradient(Vj)/(100.*dt)+10**(-30))))/(100.*dt)
	nu=-1./dlndVj
	nu2=nu[where((nu>0.) & (nu<1.))]
	maxNu.append(nu2[argmin(abs(gradient(nu2)))])
	#plot(nu2,marker='o')


fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(maxEtaZ,label=r'$\eta_Z $',marker='s',color='b')
ax1.plot(maxEtaX,label=r'$\eta _X $',marker='o',color='b')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel(r'$\eta$', color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')
ax1.legend(loc='upper center')



ax2.plot(maxNu,label=r'$\nu $',marker='x',color='r')
# Make the y-axis label and tick labels match the line color.
ax2.set_ylabel(r'$\nu$', color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')
ax2.legend(loc='upper right')



plt.show()
