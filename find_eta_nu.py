from pylab import *
from numpy import *
import  matplotlib.pyplot as plt

from load_results import *

maxEtaZ=[]
maxEtaX=[]
maxNuPos=[]
maxNuNeg=[]

NdichoPlot = matrixKappa.shape[0]

for i in range(NdichoPlot):
	EtaZt=matrixEtaZ[i]
	EtaXt=matrixEtaX[i]
	#dtEtaZ = abs(gradient(EtaZt)/dt)
	#ind = argmin(dtEtaZ)
	#maxEtaZ.append(EtaZt[ind])
	maxEtaZ.append(EtaZt[argmin(abs(gradient(EtaZt)))])
	maxEtaX.append(EtaXt[argmin(abs(gradient(EtaXt)))])


#indexVj=0
#matrixV = matrixy#[istep][:][0:Nrho]
#for istep in range(Ndicho):
#	Vj=[]
#  	NTi=len(matrixEtaZ[istep])  
##  	ti=linspace(0,100*dt*NTi,NTi)
#	ti=matrixT[istep]
#	dti=gradient(ti)
#	for k in range(len(matrixV[istep])):
#		Vj.append(matrixV[istep][k][indexVj].real)
#	dlndVj=gradient(log(abs(gradient(Vj)/(100.*dti)+10**(-30))))/(100.*dti)
#	nu=-1./dlndVj
#	nu2=nu[where((nu>0.) & (nu<1.))]
#	maxNu.append(nu2[argmin(abs(gradient(nu2)))])
#	#plot(nu2,marker='o')

indexVj=0
for istep in range(NdichoPlot):
	Vj=[]
  	ti=matrixT[istep]
	if len(ti)<2:
		break
	dti=gradient(ti)
	for k in range(len(matrixy[istep])):
		Vj.append(matrixy[istep][k][indexVj].real)
	dlndVj=gradient(log(abs(gradient(Vj,dti)+10**(-30))),dti)
	nu=-1./dlndVj
	nu2=nu[where((nu>0.) & (nu<1.))]
	dnu2=gradient(nu2)
	d2nu2=gradient(dnu2)
	positionMin = argmin(abs(dnu2))
	if d2nu2[positionMin]>0.:
		maxNuPos.append(nu2[positionMin])
	else:
		maxNuNeg.append(nu2[positionMin])
		
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



ax2.plot(maxNuPos,label=r'$\nu $',marker='x',color='r')
ax2.plot(maxNuNeg,marker='x',color='g')
# Make the y-axis label and tick labels match the line color.
ax2.set_ylabel(r'$\nu$', color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')
ax2.legend(loc='upper right')

mZ= maxEtaZ[-1]

print maxEtaZ[-1],maxEtaX[-1],(maxNuPos[-1]+maxNuNeg[-1])/2.
print 'etaZ={0:.6f}, etaX={1:.6f}, nu={2:.6f}'.format(maxEtaZ[-1],maxEtaX[-1],(maxNuPos[-1]+maxNuNeg[-1])/2.)
plt.show()
