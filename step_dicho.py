from pylab import *
import numpy as np

from global_variables import *
from time_step import *

def step_dicho(Vinit):
  
  V=Vinit
  etaZ=0.01
  etaX=0.01
  etaZPlot=[]
  etaXPlot=[]
  Vplot=[]

  for n in range(NT):
	  V,etaZ,etaX,VZXerror=time_step(V,etaZ,etaX)
	  Vp=d_rho(V)
	  if min(Vp)>0:
	    Vpre=V
	    Vprepre=Vpre
	  else:
	    print "BREAK because of non-monotonity of V' in simu beta=",beta
	    if Vpre[0]>0.:
	    	phase=1
	    break
	  if Vpre[0]<-0.7:
	  	phase=0
	  	break
	  if Vpre[0]>0:
	  	phase=1
		break

	  if(n%(100)==0):
		  etaZPlot.append(etaZ)
		  etaXPlot.append(etaX)
		  Vplot.append(V)
		if max(abs(VZXerror))>0.:
			print "V,Z or X is not real.... MAIS LOL QUOI !! VZXerror=",VZXerror
		  
	
  return phase,etaZPlot,etaXPlot,Vplot
