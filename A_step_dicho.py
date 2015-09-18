from pylab import *
import numpy as np

from global_variables import *
from A_time_step import *

def step_dicho(Vinit,Zinit,Xinit):
  
	VZX=[Vinit,Zinit,Xinit]
	etaZ=0.01
	etaX=0.01
	etaZPlot=[]
	etaXPlot=[]
	VZXplot=[]

	for n in range(NT):
		VZX,etaZ,etaX,VZXerror=time_step(VZX,etaZ,etaX)
		V=VZX[0]
		Vp=d_rho(V)		
		if min(Vp)<0:
	    		print "BREAK because of non-monotonity of V' in simu beta=",beta
	    	if V[0]>0.:
	    		phase=1
	    		break
	  	if V[0]<-0.7:
	  		phase=0
		  	break

		if(n%(100)==0):
			etaZPlot.append(etaZ)
			etaXPlot.append(etaX)
			VZXplot.append(VZX)
		if max(abs(VZXerror))>0.:
			print "V,Z or X is not real.... MAIS LOL QUOI !! VZXerror=",VZXerror
		  
	
	return phase,etaZPlot,etaXPlot,VZXplot
