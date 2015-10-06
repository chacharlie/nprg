from pylab import *
from numpy import *
from scipy import interpolate

from global_variables import beta,omegcol,Nomeg,NQ,omegaRangeDico


############# Importing R1 (from data created using Mathematica)
sbeta=str(beta)
extension='..dat'
if sbeta=='1.0':
	sbeta='1'
	extension='.dat'
elif sbeta=='5.0':
	sbeta='5'
	extension='.dat'
elif sbeta=='10.0':
	sbeta='10'
	extension='.dat'
omegRange=omegaRangeDico[sbeta]

folderPath='data/'

R1reFile=folderPath+'R1FTre_beta-'+sbeta+'_omegRange-'+str(omegRange)+extension
R1imFile=folderPath+'R1FTim_beta-'+sbeta+'_omegRange-'+str(omegRange)+extension


dataRe = loadtxt(R1reFile)
dataIm = loadtxt(R1imFile)

sizeRe=dataRe.shape[0]
sizeIm=dataIm.shape[0]

omegaRe=zeros(sizeRe)
omegaIm=zeros(sizeIm)
R1re=zeros(sizeRe)
R1im=zeros(sizeIm)

for i in range(sizeRe):
	omegaRe[i]=dataRe[i][0]
	R1re[i]=dataRe[i][1]

for i in range(sizeIm):
	omegaIm[i]=dataIm[i][0]
	R1im[i]=dataIm[i][1]


######### Interpolating R1
iR1re = interpolate.splrep(omegaRe,R1re,s=0)
iR1im = interpolate.splrep(omegaIm,R1im,s=0)

omega = omegcol.real # ATTE?TION : pour plusomega, signe -

R1col = interpolate.splev(omega,iR1re,der=0) + 1j*interpolate.splev(omega,iR1im,der=0)
domegR1col = interpolate.splev(omega,iR1re,der=1) + 1j*interpolate.splev(omega,iR1im,der=1)

R1omeg = np.zeros((Nomeg,NQ),dtype=complex)
domegR1omeg = np.zeros((Nomeg,NQ),dtype=complex)
for itemp in range(NQ):
	R1omeg[:,itemp] = R1col[:]
	domegR1omeg[:,itemp] = domegR1col[:]

