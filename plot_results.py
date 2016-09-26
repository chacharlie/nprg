from pylab import *
from numpy import *
from scipy import interpolate
import  matplotlib.pyplot as plt

from load_results import *
from decimal import Decimal

NdichoPlot = matrixKappa.shape[0]
indexPlotV = NdichoPlot-1

if printKappa:
	for i,m in enumerate(matrixKappa):
		print i,Decimal(m[0]),Decimal(m[1])

tmin=0.
for i in range(NdichoPlot):
	ti=matrixT[i]  
  	if ti[-1]<tmin:
    		tmin=ti[-1]


if model=='A' and afficheVZX:
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
#	fig,ax = plt.subplots()
	indexVj=10
	
	for istep in range(NdichoPlot):
		Vj=[]
	  	ti=np.array(matrixT[istep])
		if len(ti)<2:
			break
		dti=gradient(ti)
		for k in range(len(matrixy[istep])):
			Vj.append(matrixy[istep][k][indexVj].real)
		Vj = np.array(Vj) # correct
		dlndVj=gradient(log(abs(gradient(Vj,dti)+10**(-25))),dti) #correct

##		tgrid=linspace(0,-ti[-1],100)
##		Vtck = interpolate.splrep(-ti,Vj)
##		Vp = interpolate.splev(tgrid,Vtck,der=1)
##		VV = interpolate.splev(tgrid,Vtck,der=0)
##		logVp = log(abs(Vp+1e-25))
##		logVptck = interpolate.splrep(tgrid,logVp)
##		dlndV = interpolate.splev(tgrid,logVptck,der=1)
#       	
#		tlen=len(ti)
#		mask = np.index_exp[2000:3200]
#		tshort=-ti[mask]
#		print tshort[0],tshort[-1],tlen
#
#		Vpol = np.polyfit(tshort,Vj[mask],deg=12)
#		Vpolp = np.polyder(Vpol)
#		logV = log(abs(np.polyval(Vpolp,tshort)+1e-25))
#		logVpol = np.polyfit(tshort,logV,10)
#		dlndVpol = np.polyder(logVpol)
#		dlndV = np.polyval(dlndVpol,tshort)

      		ax6.plot(ti,-1./dlndVj,'o') #correct
#		ax6.plot(tshort,1./dlndV,'-')


	ax6.set_title('Nu')
	ax6.set_ylim([0,1])
	ax6.set_xlim([0,tmin])

#	ax.set_xlim([0,log(-tmin)])


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
