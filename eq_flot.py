from pylab import *
from scipy.interpolate import interp1d

from global_variables import *
from diff_op import *

def eqFlot(V,etaZ,etaX):
	Vdyn = zeros((V.size))
	Vp=d_rho(V)
	Vpp=d2_rho(V)
  
	s1 = -(etaZ*qR1+qdqR1+(2.-etaZ+etaX)*qdomegR1)
	s1R = (-1.+R2)*s1
	s2 = -(etaX*R2+qdqR2+(2-etaZ+etaX)*omegdomegR2)
	
	for k in range(Nrho):
		hLo= homeg+V[k]+2.*rho[k]*Vp[k]
		hTo= homeg+V[k]
		f= 3.*Vp[k]+2.*rho[k]*Vpp[k]
		
		Vdyn[k]=eqV(Vp[k],f,s1R,s2,hLo,hTo)

	Z0dyn,X0dyn = eq0XZ(V,Vp,etaZ)
		
 	Vdim=(-2.+etaZ)*V+(-2.+dim+etaZ)*rho*Vp 
	Z0dim= etaZ
	X0dim= etaX

	return Vdim+Vdyn,etaZ-(Z0dim+Z0dyn),etaX-(X0dim+X0dyn)


def eqV(Vpk,f,s1R,s2,hLo,hTo):
  num=(1.-NN)*hLo[::-1,:]**2.*hLo**2*(hTo*\
    s1R[::-1,:]+hTo[::-1,:]*(s1R-hTo*s2))*Vpk+\
    hLo*hTo[::-1,:]**2*hTo**2*s1R[::-1,:]*f+\
    hLo[::-1,:]*hTo[::-1,:]**2*hTo**2*(s1R-hLo*s2)*f
  den=hLo[::-1,:]**2*hLo**2*hTo[::-1,:]**2*hTo**2
  
  VdynTemp=4.*qdim*num/den
  VdynNew= dot(dot(VdynTemp,wQ),wOmeg)*1./(2.*pi)
  
  return VdynNew

def eq0XZ(V,Vp,etaZ):
	rho0,Vrho0,Vprho0 = findMinU(V,Vp)
	
	dtR=-qq[0,:]*(etaZ*regu[0,:]+2*qq[0,:]*regup[0,:])	
	hL= qq[0,:]*(regu[0,:]+1.)+Vrho0+2.*rho0*Vprho0
	hT= qq[0,:]*(regu[0,:]+1.)+Vrho0
	
	Z0dyn = eq0Z(Vprho0,dtR,hL,hT,rho0,hp,hpp)
	X0dyn = eq0X(Vprho0,dtR,hL,hT,rho0)

	return Z0dyn,X0dyn		

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

def findMinU(V,Vp):
	if int(abs(sum(sign(V))))!=Nrho:
		VinterpRho0=interp1d(V.real,rho)
		rho0=float(VinterpRho0(0.))
		Vinterp=interp1d(rho,V)
		Vpinterp=interp1d(rho,Vp)
		return rho0,Vinterp(rho0),Vpinterp(rho0)
	else:
		print "V=",V
		return 0,V[0],Vp[0]
	if max(abs(V.imag))>10**(-12.):
		print "ERROR: V is not real: Vmaxi=",max(abs(V.imag))
		return 0,V[0],Vp[0]


