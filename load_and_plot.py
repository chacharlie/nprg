from pylab import *
from numpy import *
import  matplotlib.pyplot as plt
#import os as os

Ndicho=1
Nomeg=200
Lomeg=50
beta=0.1
vmid=-0.1209

#currentPath=os.path.dirname(os.path.realpath(__file__))

folderPath='results/N1d3alpha2NT40000Nrho30NQ50/'

fileName=folderPath+'Veta-'+str(Ndicho)+'-'+str(Nomeg)+'-'+str(Lomeg)+'-'+str(beta)+'-'+str(vmid)

data=load(fileName+'.npz')
matrixZ=data['etaZResults']
matrixX=data['etaXResults']
matrixV=data['Vresults']

f1=plt.figure()
f2=plt.figure()
f3=plt.figure()
ax1= f1.add_subplot(111)
ax2 = f2.add_subplot(111)
ax3 = f3.add_subplot(111)
for i in range(Ndicho):
  ax1.plot(matrixZ[i], label="step="+str(i),marker='o')
  ax2.plot(matrixX[i], label="step="+str(i),marker='o')

index=0
for j in range(len(matrixV[index])):
  if (j%5)==0:
    ax3.plot(matrixV[index][j],label="time="+str(j),marker='o')

ax1.set_title('Eta Z')
ax1.set_ylim([0, 0.1])

ax2.set_title('Eta X')
ax2.set_ylim([0, 0.1])

ax1.legend(loc=4)
ax2.legend(loc=4)
ax3.legend(loc=4)
plt.show()