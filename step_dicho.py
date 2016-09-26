from pylab import *
import numpy as np

from global_variables import *
if RKadaptatif:
	from time_stepper2 import *
else:
	from simple_time_step import *


def step_dicho(yinit):
  
	y=yinit
	dty=zeros((size(y)))
	etaZ=0.01
	etaX=0.01
	if RKadaptatif:
		h=dt0

	stepCount=0
	t=0.	
	etaZPlot=[]
	etaXPlot=[]
	yplot=[]
	tPlot=[]

	while t>T: # (RG-time is negative)
		if RKadaptatif:
			h,hdid,y,dty,etaZ,etaX,rho0 = stepper(h,y,dty)
			t += hdid
		else:
			y,dty,etaZ,etaX,rho0 = simple_stepper(dt,y,dty)
			t += dt
		
		if model=='ON':
			yp=d_rho(y)
		elif model=='A':
			yp=d_rho(y[0:Nrho])
	    	if y[0]>0.:
			print 'y>0'
	    		phase=1
	    		break
	  	if y[0]<-0.7 and NN==1 and dim==3:
			print 'y<-0.7'
	  		phase=0
		  	break
	  	if y[0]<-2. and NN==1 and dim<3:
	  		phase=0
		  	break
		if y[0]<-1. and NN>1:
			phase=0
			break
		if etaZ<0.:
			phase=1
			print 'break because etaZ is negative (etaZ='+str(etaZ)+') phase set to 1 (high T)'
			break

		if etaX<0.:
			phase=1
			print 'break because etaX is negative (etaX='+str(etaZ)+') phase set to 1 (high T)'
			break

		if rho0==0.:
			print 'RHOOOOOOO=',rho0
	    		if y[0]>0.:
	    			phase=1
	  		else:
	  			phase=0
		  	break

		if model=='A' and approx=='4':
			maxZ=max(y[Nrho:2*Nrho])
			maxX=max(y[2*Nrho:3*Nrho])
			if maxZ<0.95:
				print 'Zmax=',maxZ
		    		if y[0]>0.:
		    			phase=1
		  		else:
		  			phase=0
			  	break
			
			if maxX<0.95:
				print 'Xmax=',maxX
		    		if y[0]>0.:
		    			phase=1
		  		else:
		  			phase=0
			  	break

		if(stepCount%(5)==0):
			if min(yp)<0:
	    			print "y is not monotonous at t=",t
				print "phase arbitrarily set to 0"
				phase=0
				break
			etaZPlot.append(etaZ)
			etaXPlot.append(etaX)
			yplot.append(y)
			tPlot.append(t)
		
		stepCount+=1

#		if max(abs(np.array(VZXerror)))>0.:
#			print "V,Z or X is not real.... MAIS LOL QUOI !! VZXerror=",VZXerror
		  
	if t<T:
		print 'time is over, phase is arbitrarily set to 0'
		phase=0

	return phase,etaZPlot,etaXPlot,yplot,tPlot
	

