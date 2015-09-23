from pylab import *
import numpy as np
from numpy import polynomial as poly

# physique
dim=3.			# dimension d'espace
NN=1. 			# dimension des spins

# numerique
NT = 40000		# nombre de pas de temps (de RG)
NQ = 50		# nombre de pas pour les impulsions
Nomeg = 200	# nombre de pas pour les frequences
Nrho= 30		# nombre de pas pour le potentiel
#rho0i = 14		# indice du potentiel pour l'evaluation de eta

# geometrique
T = -30.				# taille du domaine selon t
LQ = 4.2		# taille du domaine selon q
Lomeg = 50		# taille du domaine selon omega
Lrho = 0.0984*NN/(2**(-1-dim)*pi**((-dim/2))/math.gamma(dim/2))	# taille du domaine selon rho
dt = T/NT			# pas de temps
drho = Lrho/Nrho	# pas de potentiel

# parametres du regulateur
choixRegu=1		# 1 : regulateur "mou", 2 : regulateur plus violent...
alpha = 2.		# parametre du regulateur en impulsions
beta = 1.		# parametre du regulateur en frequences

# variables globales
rho = linspace(0,Lrho,Nrho)

# integration sur y
(qtemp,wtemp) = poly.legendre.leggauss(NQ)
qlign = LQ/2.*(1.+qtemp)+0*1j
q = np.zeros((Nomeg,NQ),dtype=complex)
for itemp in range(Nomeg):
	q[itemp,:] = qlign[:]
qq=q**2
qdim=q**(dim-1.)
wQ = LQ/2.*wtemp

# integration sur omega
(omegtemp,wOmegtemp) = poly.legendre.leggauss(Nomeg)
omegcol=Lomeg*omegtemp+0*1j
omeg= np.zeros((Nomeg,NQ),dtype=complex)
for itemp in range(NQ):
	omeg[:,itemp] = omegcol[:]
wOmeg=Lomeg*wOmegtemp


## regulateurs
#alpha = 2.		# parametre du regulateur en impulsions
#beta = 0.4		# parametre du regulateur en frequences
#regu = alpha/(exp(q**2)-1.)
#regup = -regu**2*exp(q**2)/alpha
#regupp = 2.*regu**3*exp(2.*q**2)/alpha**2+regup
#
## DANS R1, on change R1(omeg) en R1(-omeg) pour etre en accord avec les conventions (contradictoires) entre la def de R1 de facon causale, et les conventions "arbitraires" de la definition de R1(omeg)  
#R1 = 1j/(1j-beta*omeg)*regu
#dqR1 = 1j/(1j-beta*omeg)*regup
#dqqR1 = 1j/(1j-beta*omeg)*regupp
#domegR1 = +1j*beta/(1j-beta*omeg)**2*regu
##R1 = 1j/(1j+beta*omeg)*regu
##dqR1 = 1j/(1j+beta*omeg)*regup
##dqqR1 = 1j/(1j+beta*omeg)*regupp
##domegR1 = -1j*beta/(1j+beta*omeg)**2*regu
#
#R2 = qq*beta/(1+beta**2*omeg**2)*regu
#R21= -1.+R2
#dqR2 = beta/(1+beta**2*omeg**2)*(regup*qq+regu)
#dqqR2 = beta/(1+beta**2*omeg**2)*(regupp*qq+2.*regup)
#domegR2 = -2*beta**3*omeg/(1+beta**2*omeg**2)**2*regu*qq
#
## pour optimiser les calculs
#dqh = 1.+R1+qq*dqR1 #1.+regu[0,:]+qq[0,:]*regup[0,:]
#dqqh = 2.*dqR1+qq*dqqR1#2.*regup[0,:]+qq[0,:]*regupp[0,:]
#
#qR1=qq*R1
#qdqR1=2.*qq*qq*dqR1
#qdomegR1=qq*omeg*domegR1
#qdomegR12=qq*omeg*domegR1[::-1,:]
#qdqR2=2.*qq*dqR2
#omegdomegR2=omeg*domegR2
#
#homeg=qq*(R1+1.)+1j*omeg
#hR1= qq*domegR1
#homegp= 1j + qq*domegR1
