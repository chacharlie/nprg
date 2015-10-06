from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from time import time
from decimal import Decimal

from global_variables import *
from step_dicho import *
from weigth_dicho import *

Ndicho=30
propDicho=0.4

kmin=0.2
kmax=6.
kappa=(kmin+kmax)/2.

folderPath='results/N'+str(int(NN))+'d'+str(int(dim))+'alpha'+str(int(alpha))+'Nrho'+str(Nrho)+'NQ'+str(NQ)+'/'
fileName=folderPath+model+'-'+str(approx)+'-'+str(Ndicho)+'-'+str(atol)+'-'+str(rtol)+'-'+str(Nomeg)+'-'+str(Lomeg)\
	+'-'+str(beta)+'-'+str(kappa)+'-'+str(choixRegu)+'-'+str(propDicho)+'-moinsomega-stepper2-diff-'+str(diffOrder)+'-'+str(edgeOrder)+'-Tmax'+str(T)
	
etaZResults=[]
etaXResults=[]
yResults=[]
kappaResults=[]
tResults=[]

countPhase0=0
countPhase1=0

tsimu = time()
for i in range(Ndicho):
  print "i=",i, " time =",(time()-tsimu),"for simu ",fileName
  if (countPhase0+countPhase1)<2:
    kappa=(kmin+kmax)/2.
  else:
    wMin,wMax = weigthDicho(etaZlowT,etaZhighT,propDicho)
    kappa = wMin*kmin+wMax*kmax	
    print "weigth is ON !! SBBLAAHH !!",wMin,wMax
  
  print "kappa=",Decimal(kappa)
  
  if model=='ON':
    Vinit=0.1*(rho-kappa) 
    phase,etaZPlot,etaXPlot,yPlot,tPlot=step_dicho(Vinit) #phase=0 : low temp, =1, high temp
  elif model=='A':
    Vinit=0.1*(rho-kappa)  
    Zinit= ones((rho.size))
    #Zinit= np.array([1.00315, 1.00867, 1.0141, 1.01942, 1.02464, 1.02973, 1.03468, 1.03947, 1.0441, 1.04855, 1.05282, 1.05689, 1.06076, 1.06443, 1.06788, 1.07112, 1.07415, 1.07697, 1.07958, 1.08198, 1.08418, 1.08618, 1.08799, 1.08962, 1.09108, 1.09237, 1.0935, 1.09448, 1.09532, 1.09603])
    Xinit = ones((rho.size))
    #Xinit= np.array([0.00350089, 0.00353313, 0.00356449, 0.0035949, 0.00362429, 0.00365258, 0.0036797, 0.00370562, 0.00373026, 0.00375359, 0.00377559, 0.00379622, 0.00381548, 0.00383337, 0.00384988, 0.00386503, 0.00387885, 0.00389138, 0.00390264, 0.00391268, 0.00392156, 0.00392932, 0.00393603, 0.00394174, 0.00394652, 0.00395041, 0.0039535, 0.00395583, 0.00395747, 0.00395847])
    yInit=np.array([Vinit,Zinit,Xinit])
    phase,etaZPlot,etaXPlot,yPlot,tPlot=step_dicho(yInit.flatten()) #phase=0 : low temp, =1, high temp
  
  if phase==0:
    kmax = kappa
    countPhase0=1
    etaZlowT=etaZPlot
  else:
    kmin=kappa
    countPhase1=1
    etaZhighT=etaZPlot
  
  etaZResults.append(np.array(etaZPlot))
  etaXResults.append(np.array(etaXPlot))
  yResults.append(np.array(yPlot))
  kappaResults.append(kappa)
  tResults.append(tPlot)


etaZResults=np.array(etaZResults)
etaXResults=np.array(etaXResults)
yResults=np.array(yResults)
kappaResults=np.array(kappaResults)
tResults=np.array(tResults)
np.savez(fileName,etaZResults=etaZResults,etaXResults=etaXResults,yResults=yResults,kappaResults=kappaResults,tResults=tResults)
print "SIMU ",fileName," IS TERMINAOUCH after ",time()-tsimu,"seconds"
