from pylab import *
import numpy as np

from global_variables import *
from time_step import *

def step_dicho(Vinit):
  
  V=Vinit
  etaZ=0.1
  etaX=0.1
  phase=0
  etaZPlot=[]
  etaXPlot=[]
  Vplot=[]

  for n in range(NT):
	  V,etaZ,etaX=time_step(V,etaZ,etaX)
	  Vp=d_rho(V)
	  if min(Vp)>0:
	    Vpre=V
	    Vprepre=Vpre
	  else:
	    print "BREAK because of non-monotonity of V' in simu beta=",beta
	    break
	  #if min(Vp)<0:
	    #print "There might be a problem at t=",n*dt,"for simu beta=",beta 
	  if abs(max(V))>10.:
		  break

	  if(n%(100)==0):
		  etaZPlot.append(etaZ)
		  etaXPlot.append(etaX)
		  Vplot.append(V)
		  
  if min(Vpre[15:]-Vprepre[15:])>0:
    phase=1
  #else:
    #temp= np.sign(np.gradient(etaZPlot))
    #temp2=temp[1:]-temp[:-1]
    #temp3=np.where(temp2==0)
    #if len(temp3[0])>1:
      #phase=1
	
  return phase,etaZPlot,etaXPlot,Vplot
