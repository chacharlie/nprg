from pylab import *

from global_variables import *
from diff_op import *
from A_eq_flot import * 

def time_step(VZX,etaZ,etaX):

	#Runge-Kunta d'ordre 4
	dtVZX,etaZ_new,etaX_new,VZXerror = eqFlot(VZX,etaZ,etaX,1)
	VZX1=dtVZX

	dtVZX,_ = eqFlot(VZX+dt/2.*VZX1,etaZ,etaX,0)
	VZX2=dtVZX

	dtVZX,_ = eqFlot(VZX+dt/2.*VZX2,etaZ,etaX,0)
	VZX3=dtVZX

	dtVZX,_ = eqFlot(VZX+dt*VZX3,etaZ,etaX,0)
	VZX4=dtVZX
	
	VZX_new = VZX + dt/6.*(VZX1+2*VZX2+2*VZX3+VZX4)

	#dtZ0,dtX0=findZX0(VZX_new[0],VZX_new[1],VZX_new[2])
	#dtZ0,dtX0=findZX0(VZX_new[0],dtVZX[1],dtVZX[2])
	
	return VZX_new,etaZ_new,etaX_new,VZXerror
