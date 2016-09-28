from pylab import *
import numpy as np
from numpy import polynomial as poly

#Pour importer R1 si choixRegu==2
omegaRangeDico={'0.01':1000,'0.02':500,'0.05':400,'0.07':142.857,'0.1':100.,'0.21':47.619,'0.25':40.,'0.5':20.,'0.75':13.3333,'1.0':10,'5.0':20,'10.0':10}

# parametres du regulateur
choixRegu=1		# 1 : regulateur "mou", 2 : regulateur plus violent...
choiceReguQ=0	# 0 : regu expo, 2+ : regu litim d'ordre 2+	
alpha = 2.	# parametre du regulateur en impulsions
beta = 0		# parametre du regulateur en frequences
exact= True
aa = 1
DD = 4.

#modele
model='A'	# 'A' or 'ON'
approx=1	# (pour model A) 1: LPA', 2: Z complet, 3: X complet, 4: X,Z complets

# physique
dim=3			# dimension d'espace
NN=1 			# dimension des spins

# numerique
NQ = 30		# nombre de pas pour les impulsions
Nomeg = 200	# nombre de pas pour les frequences
Nrho= 2*30	# nombre de pas pour le potentiel

# Ordre des derivees en rho, et ordre des derivees au bord...
diffOrder = 5 
edgeOrder = 5 

T = -30.		# taille maximale du domaine selon le temps de RG

RKadaptatif=False
if RKadaptatif:
	#Runge-Kutta adaptatif
	dt0 = -1e-4		# nombre de pas de temps initial
	atol = 1.e-6	# tolerance absolue sur l'erreur dans le Runge-Kutta
	rtol = 0*1.e-6  	# tolerance relative
else:
	#pas adaptatif
	NT = 40000
	dt = T/NT

# geometrique
if choiceReguQ == 0:
	LQ = 4.2		# taille du domaine selon q
elif choiceReguQ >1:
	LQ = 1.

if choixRegu==1:
	Lomeg = 50.		# taille du domaine selon omega
elif choixRegu==2:
	Lomeg=omegaRangeDico[str(beta)]
#Lrho = 0.0984*NN/(2**(-1-dim)*pi**((-dim/2))/math.gamma(dim/2))	# taille du domaine selon rho
Lrho = 15.#1.*0.0984*NN/(2**(-1-dim)*pi**((-dim/2))/math.gamma(dim/2))	# taille du domaine selon rho
drho = Lrho/(Nrho-1)	# pas de potentiel


# variables globales
rho = linspace(0,Lrho,Nrho)

## # integration sur y
#q = np.linspace(1e-6,LQ,NQ)
#qq=q**2
#qdim=q**(dim-1.)
#wQ = LQ/float(NQ)*np.ones(NQ)


# integration sur y
(qtemp,wtemp) = poly.legendre.leggauss(NQ)
qlign = LQ/2.*(1.+qtemp)+0*1j
q = np.zeros((Nomeg,NQ),dtype=complex)
for itemp in range(Nomeg):
	q[itemp,:] = qlign[:]
if beta==0 and exact==True:
	q = LQ/2.*(1.+qtemp)

qq=q**2
qdim=q**(dim-1.)	#qdim est pratique car on a un facteur y**(dim/2)=q**dim
			#mais on a aussi facteur 2*q provenant du chgmt de variable
			#sur y, etc.
wQ = LQ/2.*wtemp

# integration sur omega
(omegtemp,wOmegtemp) = poly.legendre.leggauss(Nomeg)
omegcol=Lomeg*omegtemp+0*1j
omeg= np.zeros((Nomeg,NQ),dtype=complex)
for itemp in range(NQ):
	omeg[:,itemp] = omegcol[:]
wOmeg=Lomeg*wOmegtemp

