from pylab import *
import numpy as np
from scipy.interpolate import interp1d
import  matplotlib.pyplot as plt
from numpy import gradient as gradient

from global_variables import *
from simple_time_step import *

if model=='ON':
	from eq_flow import *
elif model=='A':
	from A_eq_flow import *

if dim==2 or dim==3 or dim==4:
	dims = str(dim)
else: 
	dims = str(int(dim*10))

Lrhos = "%.2f" % Lrho

Ndicho=30
propDicho=0.5

kmin=1.#4.#1.
kmax=6.#9.#6.
kappa=(kmin+kmax)/2.

folderPath='results/March16/N'+str(NN)+'d'+dims+'Nrho'+str(Nrho)+'NQ'+str(NQ)+'/'

str1=model+'-'+str(approx)+'-'+str(Ndicho)+'-'
str3=str(Nomeg)+'-'+str(Lomeg)+'-'+str(beta)+'-'+str(alpha)+'-'+str(kappa)+'-'+str(choixRegu)+'-'+str(propDicho)+'-'+str(diffOrder)+'-'+str(edgeOrder)+'-Tmax'+str(T)+'-'+str(choiceReguQ)+'-'+Lrhos

if RKadaptatif:
	str2=str(atol)+'-'+str(rtol)+'-'
else:	
	str2=str(NT)+'-'

filePath=str1+str2+str3

fileName = folderPath + filePath# +'simpleLrho'#+'simpleInteg' #+ 'doubledLrho'# +'modifiedV0' #+'doubledLrho'# + 'test'

data=load(fileName+'.npz')
matrixEtaZ=data['etaZResults']
matrixEtaX=data['etaXResults']
matrixy=data['yResults']

NdichoPlot = matrixEtaZ.shape[0]
indexSol = NdichoPlot-1

soly = matrixy[indexSol]
solEtaX = matrixEtaX[indexSol]
solEtaZ = matrixEtaZ[indexSol]

f1=plt.figure()
f2=plt.figure()
f3=plt.figure()
ax1 = f1.add_subplot(111)
ax2 = f2.add_subplot(111)
ax3 = f3.add_subplot(111)
ax1.set_title('V(j)-V(j-1)')
ax2.set_title('dt*eqdtV')
ax3.set_title('EtaZ')

#jdeb = int(-0.1/(5*dt))
jdeb = 1
#jfin = int(-10./(5*dt))
jfin = len(soly)
print len(soly),jdeb,jfin


for j in range(jdeb,jfin,500):
	VZX = soly[j]#-soly[j-1]
	dVZX = soly[j]-soly[j-1]
	etaX = solEtaX[j]
	etaZ = solEtaZ[j]
	dtVZX = eqFlow(VZX,etaZ,etaX,0)
	ax1.plot(rho,dVZX[:1*Nrho],label="time="+str(j),marker='o')
	ax2.plot(rho,dtVZX[:1*Nrho],label="time="+str(j),marker='o')
#	ax3.plot(j,etaZ,label="time="+str(j),marker='o')
		
	if j==jdeb:
		testQ = dot(qq*regu,wQ)
		x = (2.+0.2*rho-0.01*rho**2)
		f = np.log(x)
		fp = (0.2-0.02*rho)/x
		fpp = -fp**2-0.02/x 
		f4 = -6.*fp**4-0.24*fp**3/(0.2-0.02*rho)-0.0012/x**2 
		testRho = d_rho(f)
		testRho2 = d2_rho(f)
		rhoo=0.
		f0 = np.log(2.+0.2*rhoo-0.01*rhoo**2)
		rhoo=0.174592
		f1 = np.log(2.+0.2*rhoo-0.01*rhoo**2)
		rhoo=0.349184
		f2 = np.log(2.+0.2*rhoo-0.01*rhoo**2)
		rhoo=0.523777
		f3 = np.log(2.+0.2*rhoo-0.01*rhoo**2)
		rhoo=0.698369
		f4 = np.log(2.+0.2*rhoo-0.01*rhoo**2)
		testos = -25./12.*f0+4*f1-3*f2+4./3.*f3-f4/4. 
#		ax3.plot(rho,testRho,label="d_rho(f)",marker='s')
#		ax3.plot(rho,fp,label="f\'",marker='o')
#		ax3.plot(rho,testRho2,label="d2_rho(f)",marker='s')
#		ax3.plot(rho,fpp,label="f\'\'",marker='o')
		ax3.plot(rho,testRho-fp,label="d_rho(f)-f\'",marker='s')
		ax3.plot(rho,testRho2-fpp,label="d2_rho(f)-f\'\'",marker='o')
#		ax3.plot(rho,drho**3*f4,label="theoretical error",marker='o')
#		ax3.plot(rho,gradient(f)/drho-fp,label="grad(f)-f\'\'",marker='o')
		print testQ, drho**3,drho,testos/drho, rho[0],rho[1],rho[2],rho[3],rho[4]
		print rho

	

ax1.legend(loc='best')
ax3.legend(loc='best')


plt.show()
