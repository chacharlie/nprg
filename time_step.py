from pylab import *

from global_variables import *
from diff_op import *
from eq_dim import *
from eq_dyn import *

def time_step(V,etaZ,etaX):

	#Runge-Kunta d'ordre 4
	Vdim,Vp,Vpp,Zdim,Xdim = eqDim(V,etaZ,etaX)
	Vdyn,Z0dyn,X0dyn = eqDyn(V,Vp,Vpp,etaZ,etaX)
	V1=Vdim+Vdyn
	#etaZ1 = -(Zdim+Z0dyn)
	#etaX1 = -(Xdim+X0dyn)
	
	Vdim,Vp,Vpp,Zdim,Xdim = eqDim(V+dt/2.*V1,etaZ,etaX)
	Vdyn,Z0dyn,X0dyn = eqDyn(V+dt/2.*V1,Vp,Vpp,etaZ,etaX)
	V2=Vdim+Vdyn
	#etaZ2 =-(Zdim+Z0dyn)
	#etaX2 =-(Xdim+X0dyn)
	
	Vdim,Vp,Vpp,Zdim,Xdim = eqDim(V+dt/2.*V2,etaZ,etaX)
	Vdyn,Z0dyn,X0dyn = eqDyn(V+dt/2.*V2,Vp,Vpp,etaZ,etaX)
	V3=Vdim+Vdyn
	#etaZ3 = -(Zdim+Z0dyn)
	#etaX3 = -(Xdim+X0dyn)
	
	Vdim,Vp,Vpp,Zdim,Xdim = eqDim(V+dt*V3,etaZ,etaX)
	Vdyn,Z0dyn,X0dyn = eqDyn(V+dt*V3,Vp,Vpp,etaZ,etaX)
	V4=Vdim+Vdyn
	#etaZ4 = -(Zdim+Z0dyn)
	#etaX4 = -(Xdim+X0dyn)
	
	V_new = V + dt/6.*(V1+2*V2+2*V3+V4)
	Vp_new = d_rho(V_new)
	Vpp_new= d2_rho(V_new)

	Vdyn,Z0dyn,X0dyn = eqDyn(V_new,Vp_new,Vpp_new,etaZ,etaX)	
	etaZ_new=etaZ-(Zdim+Z0dyn)
	etaX_new=etaX-(Xdim+X0dyn)
	
	return V_new, Vp_new,Vpp_new, etaZ_new, etaX_new
