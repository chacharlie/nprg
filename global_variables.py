from pylab import *
import numpy as np
from numpy import polynomial as poly

#modele
model='A'
approx=3	#1: LPA', 2: Z complet, 3: X complet, 4: X,Z complets


# physique
dim=3.			# dimension d'espace
NN=1. 			# dimension des spins

# numerique
NQ = 50		# nombre de pas pour les impulsions
Nomeg = 400	# nombre de pas pour les frequences
Nrho= 60		# nombre de pas pour le potentiel

# Ordre des derivees en rho, et ordre des derivees au bord...
diffOrder = 5 
edgeOrder = 5 

#Runge-Kutta adaptatif
dt0 = -1e-4		# nombre de pas de temps initial
T = -30.		# taille maximale du domaine selon le temps de RG
atol = 1.e-6	# tolerance absolue sur l'erreur dans le Runge-Kutta
rtol = 0*1.e-6  	# tolerance relative


# geometrique
LQ = 4.2		# taille du domaine selon q
Lomeg = 75.		# taille du domaine selon omega
Lrho = 0.0984*NN/(2**(-1-dim)*pi**((-dim/2))/math.gamma(dim/2))	# taille du domaine selon rho
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

