from pylab import *
from numpy import *
import  matplotlib.pyplot as plt
#import os as os

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
beta=1.
kappa=3.	

#currentPath=os.path.dirname(os.path.realpath(__file__))

folderPath='results/N1d3alpha2NT40000Nrho30NQ50/'

fileName=folderPath+'Veta-'+str(Ndicho)+'-'+str(Nomeg)+'-'+str(Lomeg)+'-'+str(beta)+'-'+str(kappa)+'-plusomega-regu2'

data=load(fileName+'.npz')
matrixZ=data['etaZResults']
matrixX=data['etaXResults']
matrixV=data['Vresults']

f1=plt.figure()
f2=plt.figure()
f3=plt.figure()
ax1= f1.add_subplot(111)
ax2 = f2.add_subplot(111)
ax3 = f3.add_subplot(111)
tmin=0.
for i in range(Ndicho):
  NTi=len(matrixZ[i])  
  ti=linspace(0,100*dt*NTi,NTi)
  if ti[-1]<tmin:
	tmin=ti[-1]
  ax1.plot(ti,matrixZ[i], label="step="+str(i),marker='o')
  ax2.plot(ti,matrixX[i], label="step="+str(i),marker='o')

index=0
for j in range(len(matrixV[index])):
  if (j%5)==0:
    ax3.plot(rho,matrixV[index][j],label="time="+str(j),marker='o')


ax1.set_title('Eta Z')
ax1.set_ylim([0, 0.1])
ax1.set_xlim([0,tmin])

ax2.set_title('Eta X')
ax2.set_ylim([0, 0.1])
ax2.set_xlim([0,tmin])

ax1.legend(loc=4)
ax2.legend(loc=4)
ax3.legend(loc=4)

f4=plt.figure()
ax4=f4.add_subplot(111)

indexVj=0
for istep in range(Ndicho):
	Vj=[]
  	NTi=len(matrixZ[istep])  
  	ti=linspace(0,100*dt*NTi,NTi)
	for k in range(len(matrixV[istep])):
		Vj.append(matrixV[istep][k][indexVj].real)
	dlndVj=gradient(log(abs(gradient(Vj)/(100.*dt)+10**(-30))))/(100.*dt)
	
	ax4.plot(ti,-1./dlndVj,marker='o')

ax4.set_title('Nu')
ax4.set_ylim([0,1])
ax4.set_xlim([0,tmin])

plt.show()
