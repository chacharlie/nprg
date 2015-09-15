from pylab import *
from scipy.interpolate import interp1d

from global_variables import *
from diff_op import *

def eqFlot(V,etaZ,etaX):
	Vdyn = zeros((V.size))
	Vp=d_rho(V)
	Vpp=d2_rho(V)
  
	s1 = -(etaZ*qR1+qdqR1+(2.-etaZ+etaX)*qdomegR1)
	s12 = -(etaZ*qR1[::-1,:]+qdqR1[::-1,:]+(2.-etaZ+etaX)*qdomegR12)#s1[::-1,:] 
	s1R = (-1.+R2)*s1
	s12R = (-1.+R2[::-1,:])*s12 #s1R[::-1,:] 
	s2 = -(etaX*R2+qdqR2+(2-etaZ+etaX)*omegdomegR2)
	
	for k in range(Nrho):
		hLo= homeg+V[k]+2.*rho[k]*Vp[k]
		hTo= homeg+V[k]
		f= 3.*Vp[k]+2.*rho[k]*Vpp[k]
		
		Vdyn[k]=eqV(Vp[k],f,s1R,s12R,s2,hLo,hTo)

	Z0dyn,X0dyn = eq0XZ(V,Vp,etaZ,s1,s12,s2)
		
 	Vdim=(-2.+etaZ)*V+(-2.+dim+etaZ)*rho*Vp 
	Z0dim= etaZ
	X0dim= etaX

	return Vdim+Vdyn,etaZ-(Z0dim+Z0dyn),etaX-(X0dim+X0dyn)


def eqV(Vpk,f,s1R,s12R,s2,hLo,hTo):
  num=(1.-NN)*hLo[::-1,:]**2.*hLo**2*(hTo*\
    s12R+hTo[::-1,:]*(s1R-hTo*s2))*Vpk+\
    hLo*hTo[::-1,:]**2*hTo**2*s12R*f+\
    hLo[::-1,:]*hTo[::-1,:]**2*hTo**2*(s1R-hLo*s2)*f
  den=hLo[::-1,:]**2*hLo**2*hTo[::-1,:]**2*hTo**2
  
  VdynTemp=4.*qdim*num/den
  VdynNew= dot(dot(VdynTemp,wQ),wOmeg)*1./(2.*pi)
  
  if abs(VdynNew.imag)>10**(-12.):
  	print "Vdyn is not REAL!!!!! LOL"
  return VdynNew.real

def eq0XZ(V,Vp,etaZ,s1,s12,s2):
	rho0,Vrho0,Vprho0 = findMinU(V,Vp)
	
	hL= qq*(R1+1.)+Vrho0+2.*rho0*Vprho0
	hT= qq*(R1+1.)+Vrho0 

	Z0dyn = eq0Z(Vprho0,rho0,hL,hT,s1,s12,s2)
	X0dyn = eq0X(Vprho0,rho0,hL,hT,s1,s12,s2)
	
#	# si on veut le flot de etaZ avec beta=0
#	hLbeta0= qq[0,:]*(regu[0,:]+1.)+Vrho0+2.*rho0*Vprho0
#	hTbeta0= qq[0,:]*(regu[0,:]+1.)+Vrho0
#	dtR=-qq[0,:]*(etaZ*regu[0,:]+2*qq[0,:]*regup[0,:])	
#	hp = 1.+regu[0,:]+qq[0,:]*regup[0,:]
#	hpp = 2.*regup[0,:]+qq[0,:]*regupp[0,:]
#	Z0dyn = eq0Zbeta0(Vprho0,dtR,hLbeta0,hTbeta0,rho0,hp,hpp)

	return Z0dyn,X0dyn

def eq0Z(Vp0,rho0,hL0,hT0,s1,s12,s2):
	temp=(4*rho0*(((-((hL0[::-1,:] - 1j*omeg)*(hT0[::-1,:] - 1j*omeg)*(hL0 + 1j*omeg)*((-1j)*hT0 + omeg)*(dim*dqR2*(hL0 + 1j*omeg)*((-1j)*hT0 + omeg)*(1j*hL0[::-1,:]*(hT0 + 1j*omeg) + (hT0 - hT0[::-1,:] + 2*1j*omeg)*omeg + hL0*(1j*hT0[::-1,:] + omeg)) - 2*qq*(-(dqqR2*(hL0 + 1j*omeg)*((-1j)*hT0 + omeg)*(1j*hL0[::-1,:]*(hT0 + 1j*omeg) + (hT0 - hT0[::-1,:] + 2*1j*omeg)*omeg + hL0*(1j*hT0[::-1,:] + omeg))) + dqqh*(hL0**2*(hT0 + hT0[::-1,:]) + hL0[::-1,:]*(hT0 + 1j*omeg)**2 - (hT0 + hT0[::-1,:])*omeg**2 + hL0*(hT0**2 + 4*1j*hT0*omeg + (2*1j*hT0[::-1,:] - omeg)*omeg))*(-1 + R2)))) + dqh*(hL0[::-1,:] - 1j*omeg)*(hT0[::-1,:] - 1j*omeg)*(hL0 + 1j*omeg)*((-1j)*hT0 + omeg)*(hL0**2*(hT0 + hT0[::-1,:]) + hL0[::-1,:]*(hT0 + 1j*omeg)**2 - (hT0 + hT0[::-1,:])*omeg**2 + hL0*(hT0**2 + 4*1j*hT0*omeg + (2*1j*hT0[::-1,:] - omeg)*omeg))*(4*dqR2*qq + dim*(-1 + R2)) + 4*dqh**2*(hL0[::-1,:]**2*(hT0 + 1j*omeg)**3*(1j*hT0[::-1,:] + omeg) + hL0**3*(1j*hL0[::-1,:] + omeg)*(hT0**2 + hT0[::-1,:]**2 + hT0*(hT0[::-1,:] + 1j*omeg) - 1j*hT0[::-1,:]*omeg - omeg**2) + hL0[::-1,:]*omeg*(hT0**3*(hT0[::-1,:] - 1j*omeg) - 2*hT0*(hT0[::-1,:] - 2*1j*omeg)*omeg**2 + hT0**2*omeg*(3*1j*hT0[::-1,:] + 4*omeg) + omeg**2*(hT0[::-1,:]**2 - 2*1j*hT0[::-1,:]*omeg - 2*omeg**2)) + omeg**2*(hT0**3*((-1j)*hT0[::-1,:] - omeg) + hT0**2*(3*hT0[::-1,:] - 4*1j*omeg)*omeg + 2*hT0*omeg**2*(1j*hT0[::-1,:] + 2*omeg) + 1j*omeg**2*(-hT0[::-1,:]**2 + 2*1j*hT0[::-1,:]*omeg + 2*omeg**2)) + hL0**2*(-3*hT0**2*(hL0[::-1,:] + hT0[::-1,:] - 2*1j*omeg)*omeg + hT0**3*(1j*hT0[::-1,:] + omeg) - 3*hT0*omeg*(hL0[::-1,:]*hT0[::-1,:] + 1j*hL0[::-1,:]*omeg + 2*omeg**2) + omeg*(-3*hL0[::-1,:]*hT0[::-1,:]**2 + 3*1j*hL0[::-1,:]*hT0[::-1,:]*omeg + 3*1j*hT0[::-1,:]**2*omeg + 3*hL0[::-1,:]*omeg**2 + 4*hT0[::-1,:]*omeg**2 - 4*1j*omeg**3)) + 1j*hL0*(omeg*(-3*hT0**2*(hT0[::-1,:] - 2*1j*omeg)*omeg - 6*hT0*omeg**3 + hT0**3*(1j*hT0[::-1,:] + omeg) + omeg**2*(3*1j*hT0[::-1,:]**2 + 4*hT0[::-1,:]*omeg - 4*1j*omeg**2)) + hL0[::-1,:]*(hT0**3*(hT0[::-1,:] - 1j*omeg) + 3*1j*hT0**2*hT0[::-1,:]*omeg - 6*hT0*hT0[::-1,:]*omeg**2 + omeg**2*(-3*hT0[::-1,:]**2 + 2*1j*hT0[::-1,:]*omeg + 2*omeg**2))))*qq*(-1 + R2))*s12)/(hL0 + 1j*omeg)**3 + ((1j*hL0[::-1,:] + omeg)*(1j*hT0[::-1,:] + omeg)*(-(dim*dqh*(hL0 + 1j*omeg)*(1j*hL0[::-1,:] + omeg)**2*((-1j)*hT0 + omeg)**2) + 4*dqh**2*(1j*hL0[::-1,:] + omeg)**2*((-1j)*hT0 + omeg)**2*qq - 1j*((-1j)*hL0 + omeg)*(dim*dqh*((-1j)*hL0 + omeg)*((-1j)*hT0 + omeg)*(1j*hT0[::-1,:] + omeg)**2 + 4*dqh**2*(hL0 + 1j*omeg)*(1j*hT0[::-1,:] + omeg)**2*qq + 2*((-1j)*hT0 + omeg)*(dqqh*(1j*hL0[::-1,:] + omeg)**2*((-1j)*hT0 + omeg) + dqqh*((-1j)*hL0 + omeg)*(1j*hT0[::-1,:] + omeg)**2)*qq))*(-1 + R2)*s12)/((-1j)*hL0 + omeg)**3 + ((1j*hL0[::-1,:] + omeg)**2*(1j*hT0[::-1,:] + omeg)**2*((dim*dqh*(hT0 + 1j*omeg)*((-1j)*hL0 + omeg)*(1j*hL0[::-1,:] + omeg) - 4*dqh**2*(1j*hL0[::-1,:] + omeg)*((-1j)*hT0 + omeg)*qq + (hL0 + 1j*omeg)*(dim*dqh*((-1j)*hT0 + omeg)*(1j*hT0[::-1,:] + omeg) - 4*dqh**2*(hT0[::-1,:] - 1j*omeg)*qq + 2*dqqh*(hT0 + 1j*omeg)*(hL0[::-1,:] + hT0[::-1,:] - 2*1j*omeg)*qq))*(-1 + R2)*s1 + (dim*dqh*(hL0[::-1,:] - 1j*omeg)*(hL0 + 1j*omeg)*((-1j)*hT0 + omeg)**2 - 4*dqh**2*(hL0[::-1,:] - 1j*omeg)*((-1j)*hT0 + omeg)**2*qq + ((-1j)*hL0 + omeg)*(dim*dqh*((-1j)*hL0 + omeg)*((-1j)*hT0 + omeg)*(1j*hT0[::-1,:] + omeg) + 4*dqh**2*(hL0 + 1j*omeg)*(1j*hT0[::-1,:] + omeg)*qq + 2*((-1j)*hT0 + omeg)*(dqqh*(1j*hL0[::-1,:] + omeg)*((-1j)*hT0 + omeg) + dqqh*((-1j)*hL0 + omeg)*(1j*hT0[::-1,:] + omeg))*qq))*s2))/((-1j)*hL0 + omeg)**3)*Vp0**2)/(dim*(hL0[::-1,:] - 1j*omeg)**3*(hT0 + 1j*omeg)**3*(1j*hT0[::-1,:] + omeg)**3)


	temp2=(4*rho0*(((-((hL0[::-1,:] - 1j*omeg)*(hT0[::-1,:] - 1j*omeg)*(hL0 + 1j*omeg)*(-1j*hT0 + omeg)*(dim*dqR2*(hL0 + 1j*omeg)*(-1j*hT0 + omeg)*(1j*hL0[::-1,:]*(hT0 + 1j*omeg) + (hT0 - hT0[::-1,:] + 2*1j*omeg)*omeg + hL0*(1j*hT0[::-1,:] + omeg)) - 2*qq*(-(dqqR2*(hL0 + 1j*omeg)*(-1j*hT0 + omeg)*(1j*hL0[::-1,:]*(hT0 + 1j*omeg) + (hT0 - hT0[::-1,:] + 2*1j*omeg)*omeg + hL0*(1j*hT0[::-1,:] + omeg))) + dqqh*(hL0**2*(hT0 + hT0[::-1,:]) + hL0[::-1,:]*(hT0 + 1j*omeg)**2 - (hT0 + hT0[::-1,:])*omeg**2 + hL0*(hT0**2 + 4*1j*hT0*omeg + (2*1j*hT0[::-1,:] - omeg)*omeg))*(-1 + R2)))) + dqh*(hL0[::-1,:] - 1j*omeg)*(hT0[::-1,:] - 1j*omeg)*(hL0 + 1j*omeg)*(-1j*hT0 + omeg)*(hL0**2*(hT0 + hT0[::-1,:]) + hL0[::-1,:]*(hT0 + 1j*omeg)**2 - (hT0 + hT0[::-1,:])*omeg**2 + hL0*(hT0**2 + 4*1j*hT0*omeg + (2*1j*hT0[::-1,:] - omeg)*omeg))*(4*dqR2*qq + dim*(-1 + R2)) + 4*dqh**2*(hL0[::-1,:]**2*(hT0 + 1j*omeg)**3*(1j*hT0[::-1,:] + omeg) + hL0**3*(1j*hL0[::-1,:] + omeg)*(hT0**2 + hT0[::-1,:]**2 + hT0*(hT0[::-1,:] + 1j*omeg) - 1j*hT0[::-1,:]*omeg - omeg**2) + hL0[::-1,:]*omeg*(hT0**3*(hT0[::-1,:] - 1j*omeg) - 2*hT0*(hT0[::-1,:] - 2*1j*omeg)*omeg**2 + hT0**2*omeg*(3*1j*hT0[::-1,:] + 4*omeg) + omeg**2*(hT0[::-1,:]**2 - 2*1j*hT0[::-1,:]*omeg - 2*omeg**2)) + omeg**2*(hT0**3*(-1j*hT0[::-1,:] - omeg) + hT0**2*(3*hT0[::-1,:] - 4*1j*omeg)*omeg + 2*hT0*omeg**2*(1j*hT0[::-1,:] + 2*omeg) + 1j*omeg**2*(-hT0[::-1,:]**2 + 2*1j*hT0[::-1,:]*omeg + 2*omeg**2)) + hL0**2*(-3*hT0**2*(hL0[::-1,:] + hT0[::-1,:] - 2*1j*omeg)*omeg + hT0**3*(1j*hT0[::-1,:] + omeg) - 3*hT0*omeg*(hL0[::-1,:]*hT0[::-1,:] + 1j*hL0[::-1,:]*omeg + 2*omeg**2) + omeg*(-3*hL0[::-1,:]*hT0[::-1,:]**2 + 3*1j*hL0[::-1,:]*hT0[::-1,:]*omeg + 3*1j*hT0[::-1,:]**2*omeg + 3*hL0[::-1,:]*omeg**2 + 4*hT0[::-1,:]*omeg**2 - 4*1j*omeg**3)) + 1j*hL0*(omeg*(-3*hT0**2*(hT0[::-1,:] - 2*1j*omeg)*omeg - 6*hT0*omeg**3 + hT0**3*(1j*hT0[::-1,:] + omeg) + omeg**2*(3*1j*hT0[::-1,:]**2 + 4*hT0[::-1,:]*omeg - 4*1j*omeg**2)) + hL0[::-1,:]*(hT0**3*(hT0[::-1,:] - 1j*omeg) + 3*1j*hT0**2*hT0[::-1,:]*omeg - 6*hT0*hT0[::-1,:]*omeg**2 + omeg**2*(-3*hT0[::-1,:]**2 + 2*1j*hT0[::-1,:]*omeg + 2*omeg**2))))*qq*(-1 + R2))*s12)/(hL0 + 1j*omeg)**3 + ((1j*hL0[::-1,:] + omeg)*(1j*hT0[::-1,:] + omeg)*(-(dim*dqh*(hL0 + 1j*omeg)*(1j*hL0[::-1,:] + omeg)**2*(-1j*hT0 + omeg)**2) + 4*dqh**2*(1j*hL0[::-1,:] + omeg)**2*(-1j*hT0 + omeg)**2*qq - 1j*(-1j*hL0 + omeg)*(dim*dqh*(-1j*hL0 + omeg)*(-1j*hT0 + omeg)*(1j*hT0[::-1,:] + omeg)**2 + 4*dqh**2*(hL0 + 1j*omeg)*(1j*hT0[::-1,:] + omeg)**2*qq + 2*(-1j*hT0 + omeg)*(dqqh*(1j*hL0[::-1,:] + omeg)**2*(-1j*hT0 + omeg) + dqqh*(-1j*hL0 + omeg)*(1j*hT0[::-1,:] + omeg)**2)*qq))*(-1 + R2)*s12)/(-1j*hL0 + omeg)**3 + ((1j*hL0[::-1,:] + omeg)**2*(1j*hT0[::-1,:] + omeg)**2*((dim*dqh*(hT0 + 1j*omeg)*(-1j*hL0 + omeg)*(1j*hL0[::-1,:] + omeg) - 4*dqh**2*(1j*hL0[::-1,:] + omeg)*(-1j*hT0 + omeg)*qq + (hL0 + 1j*omeg)*(dim*dqh*(-1j*hT0 + omeg)*(1j*hT0[::-1,:] + omeg) - 4*dqh**2*(hT0[::-1,:] - 1j*omeg)*qq + 2*dqqh*(hT0 + 1j*omeg)*(hL0[::-1,:] + hT0[::-1,:] - 2*1j*omeg)*qq))*(-1 + R2)*s1 + (dim*dqh*(hL0[::-1,:] - 1j*omeg)*(hL0 + 1j*omeg)*(-1j*hT0 + omeg)**2 - 4*dqh**2*(hL0[::-1,:] - 1j*omeg)*(-1j*hT0 + omeg)**2*qq + (-1j*hL0 + omeg)*(dim*dqh*(-1j*hL0 + omeg)*(-1j*hT0 + omeg)*(1j*hT0[::-1,:] + omeg) + 4*dqh**2*(hL0 + 1j*omeg)*(1j*hT0[::-1,:] + omeg)*qq + 2*(-1j*hT0 + omeg)*(dqqh*(1j*hL0[::-1,:] + omeg)*(-1j*hT0 + omeg) + dqqh*(-1j*hL0 + omeg)*(1j*hT0[::-1,:] + omeg))*qq))*s2))/(-1j*hL0 + omeg)**3)*Vp0**2)/(dim*(hL0[::-1,:] - 1j*omeg)**3*(hT0 + 1j*omeg)**3*(1j*hT0[::-1,:] + omeg)**3)


	
	ZdynTemp=4.*qdim*temp2
	ZdynNew= dot(dot(ZdynTemp,wQ),wOmeg)*1./(2.*pi)
	
	if abs(ZdynNew.imag)>10**(-12.):
  		print "Xdyn is not REAL!!!!! LOL"
	return ZdynNew.real
	
def eq0Zbeta0(Vpk,dtR,hL0,hT0,rho0,hp0,hpp0):
	ZdynTemp=16.*qdim[0,:]/(dim*hL0**3*hT0**3)*\
	rho0*dtR*Vpk**2*(2.*hp0**2*q[0,:]**2*hT0+hL0*\
	(2.*hp0**2*qq[0,:]-hT0*(dim*hp0+2.*hpp0*qq[0,:])))
	
	return dot(ZdynTemp,wQ)
	
	
def eq0X(Vp0,rho0,hL0,hT0,s1,s12,s2):
	
	temp=(4.*rho0*((1 - 1j*hR1)*(hL0*(hT0[::-1,:] - 1j*omeg)**2 + hL0[::-1,:]**2*(hT0 + 1j*omeg) + 2*hL0[::-1,:]*omeg*(-1j*hT0 + omeg) + omeg*(1j*hT0[::-1,:]**2 + 2*hT0[::-1,:]*omeg - (hT0 + 2*1j*omeg)*omeg))*(-1 + R2)*s12 + (hL0[::-1,:]*(-1j*hT0 + omeg)**2*(domegR2*(-1j*hL0 + omeg) + 1j*(1j + hR1)*(-1 + R2)) + hL0*(2*hT0[::-1,:]*omeg*(-(domegR2*(hT0 + 1j*omeg)) + 1j*(-1 + R2) + hR1*(-1 + R2)) - omeg**2*(5 + 3*domegR2*omeg + 2*1j*hR1*(-1 + R2) + 3*1j*hR1[::-1,:]*(-1 + R2) - 5*R2) + hT0**2*(1 + domegR2*omeg + 1j*hR1[::-1,:]*(-1 + R2) - R2) + 4*hT0*omeg*(hR1[::-1,:] + 1j*(1 + domegR2*omeg - R2) - hR1[::-1,:]*R2)) + omeg*(2*omeg**2*(-1j*(2 + domegR2*omeg - 2*R2) + hR1*(-1 + R2) + hR1[::-1,:]*(-1 + R2)) - hT0*omeg*(5 + 3*domegR2*omeg + 2*1j*hR1*(-1 + R2) + 3*1j*hR1[::-1,:]*(-1 + R2) - 5*R2) + hT0[::-1,:]*omeg*(1 + domegR2*(-1j*hT0 + omeg) + 1j*hR1*(-1 + R2) - R2) + hT0**2*(2*1j + hR1 + hR1[::-1,:] + 1j*domegR2*omeg - (2*1j + hR1 + hR1[::-1,:])*R2)) + hL0**2*(hT0*(1 + domegR2*omeg + 1j*hR1[::-1,:]*(-1 + R2) - R2) + hT0[::-1,:]*(-1 + 1j*domegR2*hT0 - domegR2*omeg - 1j*hR1*(-1 + R2) + R2) + omeg*(2*1j + hR1 + hR1[::-1,:] + 1j*domegR2*omeg - (2*1j + hR1 + hR1[::-1,:])*R2)))*s12 - 1j*(1j + hR1)*(1j*hL0[::-1,:] + omeg)*(1j*hT0[::-1,:]+ omeg)*(hT0[::-1,:]*(s1 - R2*s1 + (hL0 + 1j*omeg)*s2) + hL0[::-1,:]*(s1 - R2*s1 + (hT0 + 1j*omeg)*s2) - 1j*omeg*(-2*(-1. + R2)*s1 + (hL0 + hT0 + 2*1j*omeg)*s2)))*Vp0**2)/((-1j*hL0 + omeg)**2*(1j*hL0[::-1,:] + omeg)**2*(-1j*hT0 + omeg)**2*(1j*hT0[::-1,:] + omeg)**2)

	XdynTemp=4.*qdim*temp
	XdynNew= dot(dot(XdynTemp,wQ),wOmeg)*1./(2.*pi)
	
	if abs(XdynNew.imag)>10**(-12.):
  		print "Xdyn is not REAL!!!!! LOL"
	return XdynNew.real

def findMinU(V,Vp):
	if int(abs(sum(sign(V))))!=Nrho:
		VinterpRho0=interp1d(V.real,rho)
		rho0=float(VinterpRho0(0.))
		Vinterp=interp1d(rho,V)
		Vpinterp=interp1d(rho,Vp)
		return rho0,Vinterp(rho0),Vpinterp(rho0)
	else:
		return 0,V[0],Vp[0]
	if max(abs(V.imag))>10**(-12.):
		print "ERROR: V is not real: Vmaxi=",max(abs(V.imag))
		return 0,V[0],Vp[0]


