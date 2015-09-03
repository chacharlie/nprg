from pylab import *

from global_variables import *
from diff_op import *

def eqDyn(V,Vp,Vpp,etaZ,etaX):
	Vdyn = zeros((V.size))
	
	#s11 = -(etaZ*qR1+qdqR1+(2.-etaZ+etaX)*qdomegR1)
	#s12 = -(etaZ*Rq1[::-1,:]+qdqR1[::-1,:]+(2.-etaZ+etaX)*qdomegR1[::-1,:])
	s1 = -(etaZ*qR1+qdqR1+(2.-etaZ+etaX)*qdomegR1)
	s1R = (-1.+R2)*s1
	s2 = -(etaX*R2+qdqR2+(2-etaZ+etaX)*omegdomegR2)

	dtR=-qq[0,:]*(etaZ*regu[0,:]+2*qq[0,:]*regup[0,:])	
	
	for k in range(Nrho):
		hLo= homeg+V[k]+2.*rho[k]*Vp[k]
		hTo= homeg+V[k]
		f= 3.*Vp[k]+2.*rho[k]*Vpp[k]
		
		Vdyn[k]=eqV(Vp[k],f,s1R,s2,hLo,hTo)
				
		if k==rho0i:
			hL= qq[0,:]*(regu[0,:]+1.)+V[k]+2.*rho[k]*Vp[k]
			hT= qq[0,:]*(regu[0,:]+1.)+V[k]
			Z0dyn = eq0Z(Vp[k],dtR,hL,hT,rho[k],hp,hpp)
			X0dyn = eq0X(Vp[k],dtR,hL,hT,rho[k])
		
	return Vdyn,Z0dyn,X0dyn


#def eqV(Vpk,f,s11,s12,s2,hLo,hTo):
  #num=(1.-NN)*hLo[::-1,:]**2.*hLo**2*(hTo*\
    #(-1.+R2)*s12+hTo[::-1,:]*((-1.+R2)*s11-hTo*s2))*Vpk+\
    #hLo*hTo[::-1,:]**2*hTo**2*(-1.+R2)*s12*f+\
    #hLo[::-1,:]*hTo[::-1,:]**2*hTo**2*((-1+R2)*s11-hLo*s2)*f
  #den=hLo[::-1,:]**2*hLo**2*hTo[::-1,:]**2*hTo**2
  
def eqV(Vpk,f,s1R,s2,hLo,hTo):
  num=(1.-NN)*hLo[::-1,:]**2.*hLo**2*(hTo*\
    s1R[::-1,:]+hTo[::-1,:]*(s1R-hTo*s2))*Vpk+\
    hLo*hTo[::-1,:]**2*hTo**2*s1R[::-1,:]*f+\
    hLo[::-1,:]*hTo[::-1,:]**2*hTo**2*(s1R-hLo*s2)*f
  den=hLo[::-1,:]**2*hLo**2*hTo[::-1,:]**2*hTo**2
  
  VdynTemp=4.*qdim*num/den
  VdynNew= dot(dot(VdynTemp,wQ),wOmeg)*1./(2.*pi)
  
  return VdynNew


def eq0Z(Vpk,dtR,hL0,hT0,rho0,hp0,hpp0):
	ZdynTemp=16.*qdim[0,:]/(dim*hL0**3*hT0**3)*\
	rho0*dtR*Vpk**2*(2.*hp0**2*q[0,:]**2*hT0+hL0*\
	(2.*hp0**2*qq[0,:]-hT0*(dim*hp0+2.*hpp0*qq[0,:])))
	
	return dot(ZdynTemp,wQ)
	
def eq0X(Vpk,dtR,hL0,hT0,rho0):
	XdynTemp= -8.*qdim[0,:]*rho0*dtR*Vpk**2*\
	(hL0**2+hT0**2+4.*hL0*hT0)/\
	(hL0**2*hT0**2*(hL0+hT0)**2)
	
	return dot(XdynTemp,wQ)
