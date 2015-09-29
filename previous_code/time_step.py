from pylab import *

from global_variables import *
from diff_op import *
from eq_flot import *

def time_step(V,etaZ,etaX):

	#Runge-Kunta d'ordre 4
	dtV,_,_,_ = eqFlot(V,etaZ,etaX)
	V1=dtV

	dtV,_,_,_ = eqFlot(V+dt/2.*V1,etaZ,etaX)
	V2=dtV

	dtV,_,_,_ = eqFlot(V+dt/2.*V2,etaZ,etaX)
	V3=dtV

	dtV,_,_,_ = eqFlot(V+dt*V3,etaZ,etaX)
	V4=dtV
	
	V_new = V + dt/6.*(V1+2*V2+2*V3+V4)

	Vdyn,dtZ,dtX,VZXerror = eqFlot(V_new,etaZ,etaX)	
	etaZ_new=dtZ
	etaX_new=dtX
	
	return V_new, etaZ_new, etaX_new,VZXerror

