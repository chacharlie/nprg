from pylab import *
import numpy as np
from numpy import polynomial as poly
from scipy.interpolate import interp1d

from global_variables import *
from loadR1 import R1omeg,domegR1omeg


# partie spatiale
regu = alpha/(exp(q**2)-1.)
regup = -regu**2*exp(q**2)/alpha
regupp = 2.*regu**3*exp(2.*q**2)/alpha**2+regup


########################## PARTIE FREQUENCE ##########################

###### R1 #######

# DANS R1, on change R1(omeg) en R1(-omeg) pour etre en accord avec les conventions (contradictoires) entre la def de R1 de facon causale, et les conventions "arbitraires" de la definition de R1(omeg)  
if choixRegu==1:
	R1 = 1j/(1j-beta*omeg)*regu
	dqR1 = 1j/(1j-beta*omeg)*regup
	dqqR1 = 1j/(1j-beta*omeg)*regupp
	domegR1 = +1j*beta/(1j-beta*omeg)**2*regu
	#R1 = 1j/(1j+beta*omeg)*regu
	#dqR1 = 1j/(1j+beta*omeg)*regup
	#dqqR1 = 1j/(1j+beta*omeg)*regupp
	#domegR1 = -1j*beta/(1j+beta*omeg)**2*regu

elif choixRegu==2:
	R1 = regu*R1omeg
	dqR1 = regup*R1omeg
	dqqR1 = regupp*R1omeg
	domegR1 = regu*domegR1omeg

###### R2 #######
if choixRegu==1:
	R2 = qq*beta/(1+beta**2*omeg**2)*regu
	R21= -1.+R2
	dqR2 = beta/(1+beta**2*omeg**2)*(regup*qq+regu)
	dqqR2 = beta/(1+beta**2*omeg**2)*(regupp*qq+2.*regup)
	domegR2 = -2*beta**3*omeg/(1+beta**2*omeg**2)**2*regu*qq

elif choixRegu==2:
	R2omeg = (R1omeg-R1omeg[::-1,:])/(2.*1j*omeg)
	R2 = qq*regu*R2omeg
	R21 =-1.+R2
	dqR2 = qq*(regup*qq+regu)*R2omeg
	dqqR2 =  qq*(regupp*qq+2.*regup)*R2omeg
	domegR2 =  qq*regu*(domegR1omeg+domegR1omeg[::-1,:]-(R1omeg-R1omeg[::-1,:])/omeg)/(2.*1j*omeg)


# pour optimiser les calculs
dqhnoZ = R1+qq*dqR1 #1.+regu[0,:]+qq[0,:]*regup[0,:]
dqqh = 2.*dqR1+qq*dqqR1#2.*regup[0,:]+qq[0,:]*regupp[0,:]

qdomegR1=qq*domegR1

qR1=qq*R1
qdqR1=2.*qq*qq*dqR1
qomegdomegR1=qq*omeg*domegR1
qdomegR12=qq*omeg*domegR1[::-1,:]
qdqR2=2.*qq*dqR2
omegdomegR2=omeg*domegR2

