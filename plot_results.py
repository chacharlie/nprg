from pylab import *
from numpy import *
import  matplotlib.pyplot as plt

from load_results import *
from decimal import Decimal

NdichoPlot = matrixKappa.shape[0]
indexPlotV = NdichoPlot-1

if printKappa:
	for i,m in enumerate(matrixKappa):
		print i,Decimal(m[0]),Decimal(m[1])

if model=='A':
	f4=plt.figure()
	f5=plt.figure()
	ax4= f4.add_subplot(111)
	ax5= f5.add_subplot(111)
	ax4.set_title('Z')
	ax5.set_title('X')


if afficheEta:
	f1=plt.figure()
	f2=plt.figure()
	ax1= f1.add_subplot(111)
	ax2 = f2.add_subplot(111)
	
	tmin=0.
	for i in range(NdichoPlot):
		ti=matrixT[i]  
	  	if ti[-1]<tmin:
	    		tmin=ti[-1]
		ax1.plot(ti,matrixEtaZ[i], label=str(i),marker='o')
		ax2.plot(ti,matrixEtaX[i], label=str(i),marker='o')
	
	ax1.set_title('Eta Z')
	#ax1.set_ylim([0, 0.1])
	ax1.set_xlim([0,tmin])
	ax2.set_title('Eta X')
	#ax2.set_ylim([0, 0.2])
	ax2.set_xlim([0,tmin])
	
	ax1.legend(loc='best')
	ax2.legend(loc='best')

if afficheVZX:
	f3=plt.figure()
	ax3 = f3.add_subplot(111)
	ax3.set_title(filePath)
	
	for j in range(len(matrixy[indexPlotV])):
		if (j%5)==0:
	    		ax3.plot(rho,matrixy[indexPlotV][j][0:Nrho],label="time="+str(j),marker='o')
			if model=='A':
	    			ax4.plot(rho,matrixy[indexPlotV][j][Nrho:2*Nrho],label="time="+str(j),marker='o')
	    			ax5.plot(rho,matrixy[indexPlotV][j][2*Nrho:3*Nrho],label="time="+str(j),marker='o')
	
	ax3.legend(loc='best')


if afficheNu:
	f6=plt.figure()
	ax6=f6.add_subplot(111)
	indexVj=10
	for istep in range(NdichoPlot):
		Vj=[]
	  	ti=matrixT[istep]
		if len(ti)<2:
			break
		dti=gradient(ti)
		for k in range(len(matrixy[istep])):
			Vj.append(matrixy[istep][k][indexVj].real)
		dlndVj=gradient(log(abs(gradient(Vj,dti)+10**(-25))),dti)
		
		ax6.plot(ti,-1./dlndVj,marker='o')
	
	ax6.set_title('Nu')
	ax6.set_ylim([0,1])
	ax6.set_xlim([0,tmin])

#if afficheNu==True:
#	f6=plt.figure()
#	ax6=f6.add_subplot(111)
#
#	indexVj=0
#	for istep in range(NdichoPlot):
#		Vj=[]
#	  	NTi=len(matrixEtaZ[istep])  
#	  	ti=linspace(0,100*dt*NTi,NTi)
#		for k in range(len(matrixV[istep])):
#			Vj.append(matrixV[istep][k][indexVj].real)
#		dlndVj=gradient(log(abs(gradient(Vj)/(100.*dt)+10**(-30))))/(100.*dt)
#		
#		ax6.plot(ti,-1./dlndVj,marker='o')
#	
#	ax6.set_title('Nu')
#	ax6.set_ylim([0,1])
#	ax6.set_xlim([0,tmin])

plt.show()
