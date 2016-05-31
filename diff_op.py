from pylab import *
from numpy import diff as diff
from numpy import gradient as gradient

from global_variables import *

def d_rho(F):
	Fp = zeros((F.__len__()))	
	#Fp[0] = -(Decimal(25.)*F[0]-Decimal(48.)*F[1]+Decimal(36.)*F[2]-Decimal(16.)*F[3]+Decimal(3.)*F[4])/Decimal(12.);
	#Fp[1] = -(Decimal(3.)*F[0]+Decimal(10.)*F[1]-Decimal(18.)*F[2]+Decimal(6.)*F[3]-F[4])/Decimal(12.);
	#Fp[-1]=	(Decimal(25.)*F[-1]-Decimal(48.)*F[-2]+Decimal(36.)*F[-3]-Decimal(16.)*F[-4]+Decimal(3.)*F[-5])/Decimal(12.)
	#Fp[-2] = (Decimal(3.)*F[-1]+Decimal(10.)*F[-2]-Decimal(18.)*F[-3]+Decimal(6.)*F[-4]-F[-5])/Decimal(12.)
	#Fp[2:-2] = Decimal(1./12)*F[0:-4]-Decimal(2./3)*F[1:-3]+Decimal(2./3)*F[3:-1]-Decimal(1./12)*F[4:]
	
	
	if diffOrder==3:
		Fp=gradient(F)
	elif diffOrder==5:
		Fp[2:-2] = 1./12*F[:-4]-2./3*F[1:-3]+2./3*F[3:-1]-1./12*F[4:]

	if edgeOrder==5:
		Fp[0] = -(25.*F[0]-48.*F[1]+36.*F[2]-16.*F[3]+3.*F[4])/12.
		Fp[1] = -(3.*F[0]+10.*F[1]-18.*F[2]+6.*F[3]-F[4])/12.
		Fp[-1]=	(25.*F[-1]-48.*F[-2]+36.*F[-3]-16.*F[-4]+3.*F[-5])/12.
		Fp[-2] = (3.*F[-1]+10.*F[-2]-18.*F[-3]+6.*F[-4]-F[-5])/12.

	return Fp/drho

	
def d2_rho(F):
	Fpp = zeros((F.__len__()))	



	if diffOrder==3:
		Fpp=gradient(gradient(F))
	elif diffOrder==5:
		Fpp[2:-2] = -F[0:-4]/12.+4./3*F[1:-3]-5./2*F[2:-2]+4./3*F[3:-1]-1./12*F[4:]

	if edgeOrder==5:
		Fpp[0] = 35./12*F[0]-26./3*F[1]+19./2*F[2]-14./3*F[3]+11./12*F[4]
		Fpp[1] = 11./12*F[0]-5./3*F[1]+F[2]/2.+F[3]/3.-F[4]/12.
		Fpp[-1] = 11./12*F[-5]-14./3*F[-4] +19./2*F[-3] -26./3*F[-2]+ 35./12*F[-1]
		Fpp[-2] = -F[-5]/12. + 1./3*F[-4]+F[-3]/2. -5./3*F[-2] + 11./12*F[-1]

	return Fpp/drho**2


def test_der(F):
	Fp = -(25.*F[0]-48.*F[1]+36.*F[2]-16.*F[3]+3.*F[4])/12.
	return Fp/drho

