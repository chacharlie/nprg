from pylab import *

from global_variables import *
from diff_op import *

def eqDim(V,etaZ,etaX):
	Vp = d_rho(V)
	Vpp = d2_rho(V)
	
	Vdim = (-2.+etaZ)*V[:] + (-2.+dim+etaZ)*rho[:]*Vp[:]
	Zdim = etaZ
	Xdim = etaX
	
	return Vdim,Vp,Vpp,Zdim,Xdim
