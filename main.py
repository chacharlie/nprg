from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from time import time

from global_variables import *
from step_dicho import *

Ndicho=1

vhigh=-0.1
vlow=-0.1
vmid=(vhigh+vlow)/2.

# Pour beta=0
#vhigh=-0.026953125
#vlow=-0.0271484375

folderPath='results/N'+str(int(NN))+'d'+str(int(dim))+'alpha'+str(int(alpha))+'NT'+str(NT)+'Nrho'+str(Nrho)+'NQ'+str(NQ)+'/'
fileName=folderPath+'Veta-'+str(Ndicho)+'-'+str(Nomeg)+'-'+str(Lomeg)+'-'+str(beta)+'-'+str(vmid)
etaZResults=[]
etaXResults=[]
Vresults=[]
vmidResults=[]

tsimu = time()
for i in range(Ndicho):
  print "i=",i, " time =",(time()-tsimu),"for simu ",fileName
  vmid=(vhigh+vlow)/2.
  print "vmid=",vmid
  Vinit=linspace(-0.3,0.7,Nrho)+vmid
  
  phase,etaZPlot,etaXPlot,Vplot=step_dicho(Vinit) #phase=0 : low temp, =1, high temp
  
  if phase==0:
    vlow = vmid
  else:
    vhigh=vmid
  
  etaZResults.append(np.array(etaZPlot))
  etaXResults.append(np.array(etaXPlot))
  Vresults.append(np.array(Vplot))
  vmidResults.append(vmid)


etaZResults=np.array(etaZResults)
etaXResults=np.array(etaXResults)
Vresults=np.array(Vresults)
vmidResults=np.array(vmidResults)
np.savez(fileName,etaZResults=etaZResults,etaXResults=etaXResults,Vresults=Vresults,vmidResults=vmidResults)
print "SIMU beta=",beta," IS TERMINAOUCH after ",time()-tsimu,"seconds"
