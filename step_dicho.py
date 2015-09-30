from pylab import *
import numpy as np

from global_variables import *
from time_stepper2 import *

def step_dicho(yinit):
  
	y=yinit
	dty=zeros((size(y)))
	etaZ=0.01
	etaX=0.01
	h=T/NT

	stepCount=0
	t=0.	
	etaZPlot=[]
	etaXPlot=[]
	yplot=[]
	tPlot=[]

	while t>T: # (RG-time is negative)
		h,hdid,y,dty,etaZ,etaX = stepper(h,y,dty,etaZ,etaX)
		t += hdid
		if model=='ON':
			yp=d_rho(y)
		elif model=='A':
			yp=d_rho(y[0:Nrho])
	    	if y[0]>0.:
	    		phase=1
	    		break
	  	if y[0]<-0.7:
	  		phase=0
		  	break

		if(stepCount%(5)==0):
			if min(yp)<0:
	    			print "y is not monotonous at t=",t
			etaZPlot.append(etaZ)
			etaXPlot.append(etaX)
			yplot.append(y)
			tPlot.append(t)
		
		stepCount+=1

#		if max(abs(np.array(VZXerror)))>0.:
#			print "V,Z or X is not real.... MAIS LOL QUOI !! VZXerror=",VZXerror
		  
	
	return phase,etaZPlot,etaXPlot,yplot,tPlot
