from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from time import time

from global_variables import *
from A_step_dicho import *

Ndicho=1

kmin=3.
kmax=3.
kappa=(kmin+kmax)/2.

folderPath='A_results/N'+str(int(NN))+'d'+str(int(dim))+'alpha'+str(int(alpha))+'NT'+str(NT)+'Nrho'+str(Nrho)+'NQ'+str(NQ)+'/'
fileName=folderPath+'Veta-'+str(Ndicho)+'-'+str(Nomeg)+'-'+str(Lomeg)+'-'+str(beta)+'-'+str(kappa)+'-'+str(choixRegu)+'-plusomega'
etaZResults=[]
etaXResults=[]
Vresults=[]
Zresults=[]
Xresults=[]
kappaResults=[]

tsimu = time()
for i in range(Ndicho):
  print "i=",i, " time =",(time()-tsimu),"for simu ",fileName
  kappa=(kmin+kmax)/2.
  print "kappa=",kappa
  Vinit= 0.1*(rho-kappa)
  Zinit= 0.*ones((rho.size))
  Xinit=0.1*Zinit
  
  phase,etaZPlot,etaXPlot,VZXplot=step_dicho(Vinit,Zinit,Xinit) #phase=0 : low temp, =1, high temp
  
  if phase==0:
    kmax = kappa
  else:
    kmin=kappa
  
  etaZResults.append(np.array(etaZPlot))
  etaXResults.append(np.array(etaXPlot))
  Vresults.append(np.array(VZXplot[0]))
  Zresults.append(np.array(VZXplot[1]))
  Xresults.append(np.array(VZXplot[2]))
  kappaResults.append(kappa)


etaZResults=np.array(etaZResults)
etaXResults=np.array(etaXResults)
Vresults=np.array(Vresults)
Zresults=np.array(Zresults)
Xresults=np.array(Xresults)
kappaResults=np.array(kappaResults)
np.savez(fileName,etaZResults=etaZResults,etaXResults=etaXResults,Vresults=Vresults,Zresults=Zresults,Xresults=Xresults,kappaResults=kappaResults)
print "SIMU beta=",beta," IS TERMINAOUCH after ",time()-tsimu,"seconds"
