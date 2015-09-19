from pylab import *
from numpy import *
import  matplotlib.pyplot as plt

# physique
dim=3.                  # dimension d'espace
NN=1.                   # dimension des spins

# numerique
NT = 40000              # nombre de pas de temps (de RG)
NQ = 50         # nombre de pas pour les impulsions
Nomeg = 200     # nombre de pas pour les frequences
Nrho= 30                # nombre de pas pour le potentiel

# geometrique
T = -30.                                # taille du domaine selon t
LQ = 4.2                # taille du domaine selon q
Lomeg = 50              # taille du domaine selon omega
Lrho = 0.0984*NN/(2**(-1-dim)*pi**((-dim/2))/math.gamma(dim/2)) # taille du domaine selon rho
dt = T/NT                       # pas de temps
drho = Lrho/Nrho        # pas de potentiel

# variables globales
rho = linspace(0,Lrho,Nrho)


Ndicho=1
beta=0.
kappa=3.	

choixRegu=1

folderPath='A_results/N1d3alpha2NT40000Nrho30NQ50/'

fileName=folderPath+'Veta-'+str(Ndicho)+'-'+str(Nomeg)+'-'+str(Lomeg)+'-'+str(beta)+'-'+str(kappa)+'-'+str(choixRegu)+'-moinsomega'

data=load(fileName+'.npz')
matrixEtaZ=data['etaZResults']
matrixEtaX=data['etaXResults']
matrixV=data['Vresults']
matrixZ=data['Zresults']
matrixX=data['Xresults']

print matrixV[0]
print matrixZ[0]


f1=plt.figure()
f2=plt.figure()
f3=plt.figure()
f4=plt.figure()
f5=plt.figure()
f6=plt.figure()
ax1= f1.add_subplot(111)
ax2 = f2.add_subplot(111)
ax3 = f3.add_subplot(111)
ax4 = f4.add_subplot(111)
ax5 = f5.add_subplot(111)
ax6 = f5.add_subplot(111)
tmin=0.
for i in range(Ndicho):
  NTi=len(matrixEtaZ[i])  
  ti=linspace(0,100*dt*NTi,NTi)
  if ti[-1]<tmin:
	tmin=ti[-1]
  ax1.plot(ti,matrixEtaZ[i], label="step="+str(i),marker='o')
  ax2.plot(ti,matrixEtaX[i], label="step="+str(i),marker='o')

index=0
for j in range(len(matrixV[index])):
  if (j%5)==0:
    ax3.plot(rho,matrixV[index][j],label="time="+str(j),marker='o')
    ax4.plot(rho,matrixZ[index][j],label="time="+str(j),marker='o')
    ax5.plot(rho,matrixX[index][j],label="time="+str(j),marker='o')


ax1.set_title('Eta Z')
ax1.set_ylim([0, 0.1])
ax1.set_xlim([0,tmin])

ax2.set_title('Eta X')
ax2.set_ylim([0, 0.1])
ax2.set_xlim([0,tmin])

ax1.legend(loc='best')
ax2.legend(loc='best')
ax3.legend(loc='best')
ax4.legend(loc='best')
ax5.legend(loc='best')

#indexVj=0
#for istep in range(Ndicho):
#	Vj=[]
#  	NTi=len(matrixEtaZ[istep])  
#  	ti=linspace(0,100*dt*NTi,NTi)
#	for k in range(len(matrixV[istep])):
#		Vj.append(matrixV[istep][k][indexVj].real)
#	dlndVj=gradient(log(abs(gradient(Vj)/(100.*dt)+10**(-30))))/(100.*dt)
#	
#	ax6.plot(ti,-1./dlndVj,marker='o')

ax6.set_title('Nu')
ax6.set_ylim([0,1])
ax6.set_xlim([0,tmin])

plt.show()
