from pylab import *

from global_variables import model

if model=='ON':
	from eq_flow import *
elif model=='A':
	from A_eq_flow import *

def simple_stepper(dt,y,dty):

	etaZ,etaX,rho0 = computeEta(y)

	#Runge-Kunta d'ordre 4
	k1 = eqFlow(y,etaZ,etaX,0)

	ytemp = y+dt/2.*k1
	k2 = eqFlow(ytemp,etaZ,etaX,0)

	ytemp = y+dt/2.*k2
	k3 = eqFlow(ytemp,etaZ,etaX,0)

	ytemp = y+dt*k3
	k4 = eqFlow(ytemp,etaZ,etaX,0)
	
	y_new = y + dt/6.*(k1+2*k2+2*k3+k4)

	dty_new = eqFlow(y_new,etaZ,etaX,0)	
	
	return y_new,dty_new,etaZ,etaX,rho0

def computeEta(y):
	# etaZ_new = A + B*etaZ + C*etaX
	_,A,Ap,rho0 = eqFlow(y,0.,0.,1)
	_,etaZ10,etaX10,_ = eqFlow(y,1.,0.,1)
	_,etaZ01,etaX01,_ = eqFlow(y,0.,1.,1)

	B,Bp = etaZ10-A, etaX10-Ap
	C,Cp = etaZ01-A, etaX01-Ap
	
	etaZ = (A*(1.-Cp)+C*A)/((1.-B)*(1.-Cp)-C*Bp)
	etaX = (Ap+Bp*etaZ)/(1.-Cp)

	return etaZ, etaX,rho0

