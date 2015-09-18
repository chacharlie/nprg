from pylab import *
from scipy.interpolate import interp1d

from global_variables import *
from diff_op import *
from regu import *

def eqFlot(V,Z,X,etaZ,etaX):
	Vdyn = zeros((rho.size))
	Verror= zeros((rho.size))
	Vp=d_rho(V)
	Vpp=d2_rho(V)

	Xdyn = zeros((rho.size))
	Xerror= zeros((rho.size))
	Xp=d_rho(X)
	Xpp=d2_rho(X)
	
	Zdyn = zeros((rho.size))
	Zerror= zeros((rho.size))
	Zp=d_rho(Z)
	Zpp=d2_rho(Z)


	s1 = -(etaZ*qR1+qdqR1+(2.-etaZ+etaX)*qomegdomegR1)
	s12 = -(etaZ*qR1[::-1,:]+qdqR1[::-1,:]+(2.-etaZ+etaX)*qdomegR12)#s1[::-1,:] 
	s1R = (-1.+R2)*s1
	s12R = (-1.+R2[::-1,:])*s12 #s1R[::-1,:] 
	s2 = -(etaX*R2+qdqR2+(2-etaZ+etaX)*omegdomegR2)
	
	for k in range(Nrho):
		h = qq*(R1+Z[k])+V[k]+2.*rho[k]*Vp[k]
		homeg = h+1j*omeg*X[k] #homeg+V[k]+2.*rho[k]*Vp[k]
		f= qq*Zp[k]+3.*Vp[k]+2.*rho[k]*Vpp[k]
		g = homeg*Xp[k]-X[k]*f
	
		Vdyn[k],Verror[k]=eqV(f,g,s1,s12,s2,homeg)
		Zdyn[k],Xdyn[k],Zerror[k],Xerror[k] =eqXZ(V,Vp,s1,s12,s2,rho[k],X[k],Xp[k],Xpp[k],Z[k],Zp[k],Zpp[k],f,h,homeg)

	Z0dyn,X0dyn,ZXerror = eq0XZ(V,Vp,etaZ,s1,s12,s2)
		
 	Vdim=(-2.+etaZ)*V+(-2.+dim+etaZ)*rho*Vp 
	Zdim= etaZ*Z+(-2.+dim+etaZ)*rho*Zp
	Xdim=  etaX*X+(-2.+dim+etaZ)*rho*Xp
	
	VZXerror = [max(abs(Verror))]+[max(abs(Xerror))]+[max(abs(Zerror))]

	#return Vdim+Vdyn,etaZ-(Z0dim+Z0dyn),etaX-(X0dim+X0dyn),VZXerror
	return Vdim+Vdyn,Zdim+Zdyn,Xdim+Xdyn,VZXerror


def eqV(f,g,s1,s12,s2,homeg):
	term1 = (g[::-1,:] + f*R2)*s1/(homeg[::-1,:]*homeg**2)
	term2 = (g + f*R2)*s12/(homeg[::-1,:]**2*homeg)
	term3 = -f*s2/(homeg[::-1,:]*homeg)
  
  
	VdynTemp=4.*qdim*(term1+term2+term3)
	VdynNew= dot(dot(VdynTemp,wQ),wOmeg)*1./(2.*pi)
  
	Verror=0.
	if abs(VdynNew.imag)>10**(-10.):
		Verror=VdynNew.imag
	return VdynNew.real,Verror

def eqXZ(V,Vp,s1,s12,s2,rhok,Xk,Xpk,Xppk,Zk,Zpk,Zppk,f,h,homeg):

	Xyn,Xerror = eqX(rhok,Xk,Xpk,Xppk,s1,s12,s2,f,h,homeg)
	Zdyn,Zerror = eqZ(rhok,Xk,Xpk,Zk,Zpk,Zppk,s1,s12,s2,f,homeg)

	rho0,Vrho0,Vprho0 = findMinU(V,Vp)
	Z0dyn,Zerror = eq0Z(Vprho0,rho0,hL,hT,s1,s12,s2)
	X0dyn,Xerror = eq0X(Vprho0,rho0,hL,hT,s1,s12,s2)
	
	ZXerror=[Zerror,Xerror]
#	# si on veut le flot de etaZ avec beta=0
#	hLbeta0= qq[0,:]*(regu[0,:]+1.)+Vrho0+2.*rho0*Vprho0
#	hTbeta0= qq[0,:]*(regu[0,:]+1.)+Vrho0
#	dtR=-qq[0,:]*(etaZ*regu[0,:]+2*qq[0,:]*regup[0,:])	
#	hp = 1.+regu[0,:]+qq[0,:]*regup[0,:]
#	hpp = 2.*regup[0,:]+qq[0,:]*regupp[0,:]
#	Z0dyn = eq0Zbeta0(Vprho0,dtR,hLbeta0,hTbeta0,rho0,hp,hpp)

	return Z0dyn,X0dyn,ZXerror

def eqZ(rhok,Xk,Xpk,Zk,Zpk,Zppk,s1,s12,s2,f,homeg):
	dqh = Zk + dqhnoZ	
	ZPP = Zpk + 2*rhok*Zppk
	fXp = f + 1j*omeg*Xpk
	s22 = s1*Xk+s2*homeg-R2*s1
	R2X= R2-Xk
	
	nume1 = 4*rhok*s12*R2X*(-(dim*homeg*(-(dqh*f*fXp) + (f + fXp)*homeg*Zpk)) + y*(-4*dqh**2*f*fXp + 2*dqh*(f + fXp)*homeg*Zpk + homeg*(2*dqqh*f*fXp - homeg*Zpk**2))) - dim*homeg**3*homeg[::-1,:]*s2*ZPP
	deno1 = homeg**4*homeg[::-1,:]**2

	nume2 = 4*homeg[::-1,:]*rhok*s1*Xpk*(dim*homeg*(dqh*fXp - homeg*Zpk) + 2*y*(-2*dqh**2*fXp + dqqh*fXp*homeg + dqh*homeg*Zpk)) - 4*rhok*s22*y*(-4*dqh**2*f*fXp + 2*dqh*(f + fXp)*homeg*Zpk + homeg*(2*dqqh*f*fXp - homeg*Zpk**2)) - dim*homeg*(4*dqh*f*fXp*rhok*s22 - homeg*(4*f*rhok*s22*Zpk + 4*fXp*rhok*s22*Zpk + homeg*s1*R2X*ZPP))
	deno2 = homeg**5*homeg[::-1,:]
	
	nume3 = -(s12*(4*rhok*y*(4*dqh[::-1,:]**2*fXp[::-1,:]*homeg**2*(f*R2X + homeg*Xpk) - 2*dqh[::-1,:]*homeg*homeg[::-1,:]*(2*dqR2*f*fXp[::-1,:]*homeg - 2*dqh*f*fXp[::-1,:]*R2X + homeg*((f + fXp[::-1,:])*R2X + homeg*Xpk)*Zpk) + homeg[::-1,:]*(4*dqh**2*f*homeg[::-1,:]*(fXp[::-1,:]*R2X + homeg[::-1,:]*Xpk) - 2*dqh*homeg*homeg[::-1,:]*(2*dqR2*f*fXp[::-1,:] + ((f + fXp[::-1,:])*R2X + homeg[::-1,:]*Xpk)*Zpk) + homeg*(2*dqqR2*f*fXp[::-1,:]*homeg*homeg[::-1,:] - 2*dqqh[::-1,:]*fXp[::-1,:]*homeg*(f*R2X + homeg*Xpk) + homeg[::-1,:]*(-2*dqqh*f*(fXp[::-1,:]*R2X + homeg[::-1,:]*Xpk) + homeg*Zpk*(2*dqR2*(f + fXp[::-1,:]) + R2X*Zpk))))) + dim*homeg*homeg[::-1,:]*(4*dqR2*f*fXp[::-1,:]*homeg*homeg[::-1,:]*rhok - 4*dqh[::-1,:]*fXp[::-1,:]*homeg*rhok*(f*R2X + homeg*Xpk) + homeg[::-1,:]*(-4*dqh*f*rhok*(fXp[::-1,:]*R2X + homeg[::-1,:]*Xpk) + 4*homeg*rhok*((f + fXp[::-1,:])*R2X + (homeg + homeg[::-1,:])*Xpk)*Zpk - homeg*homeg[::-1,:]*R2X*ZPP))))
	deno3 = homeg[::-1,:]**5*homeg
	
	ZdynTemp=4.*qdim*dim*(nume1/deno1+nume2/deno2+nume3/deno3)
	ZdynNew= dot(dot(ZdynTemp,wQ),wOmeg)*1./(2.*pi)
	
	Zerror=0.
	if abs(ZdynNew.imag)>10**(-10.):
		Zerror=ZdynNew.imag
	return ZdynNew.real,Zerror
	
def eqX(rhok,Xk,Xpk,Xppk,s1,s12,s2,f,h,homeg):
	XPP = Xpk + 2*rhok*Xppk
	domegh=qdomegR1

	num1 = s12*(4.*1j*domegh*f**2*h[::-1,:]*R2*rhok - 4.*1j*domegh*f**2*h[::-1,:]*rhok*Xk + 4*f**2*h*R2*rhok*Xk - 4*f**2*h[::-1,:]*R2*rhok*Xk + 4*domegh*f**2*omeg*R2*rhok*Xk - 4*f**2*h*rhok*Xk**2 + 4*f**2*h[::-1,:]*rhok*Xk**2 - 4*domegh*f**2*omeg*rhok*Xk**2 + 8.*1j*f**2*omeg*R2*rhok*Xk**2 - 8.*1j*f**2*omeg*rhok*Xk**3 - h*h[::-1,:]**2*R2*Xpk + 4.*1j*domegh*f*h[::-1,:]**2*rhok*Xpk + 4*f*h*h[::-1,:]*R2*rhok*Xpk + 4*domegh*f*h[::-1,:]*omeg*R2*rhok*Xpk + h*h[::-1,:]**2*Xk*Xpk + 2.*1j*h*h[::-1,:]*omeg*R2*Xk*Xpk - 1j*h[::-1,:]**2*omeg*R2*Xk*Xpk + 4*f*h**2*rhok*Xk*Xpk - 4*f*h*h[::-1,:]*rhok*Xk*Xpk - 4*f*h[::-1,:]**2*rhok*Xk*Xpk + 4*domegh*f*h[::-1,:]*omeg*rhok*Xk*Xpk - 8.*1j*f*h*omeg*R2*rhok*Xk*Xpk + 8.*1j*f*h[::-1,:]*omeg*R2*rhok*Xk*Xpk - 4.*1j*domegh*f*omeg**2*R2*rhok*Xk*Xpk - 2.*1j*h*h[::-1,:]*omeg*Xk**2*Xpk + 1j*h[::-1,:]**2*omeg*Xk**2*Xpk + h*omeg**2*R2*Xk**2*Xpk - 2*h[::-1,:]*omeg**2*R2*Xk**2*Xpk + 16.*1j*f*h*omeg*rhok*Xk**2*Xpk + 12*f*omeg**2*R2*rhok*Xk**2*Xpk - h*omeg**2*Xk**3*Xpk + 2*h[::-1,:]*omeg**2*Xk**3*Xpk + 1j*omeg**3*R2*Xk**3*Xpk - 12*f*omeg**2*rhok*Xk**3*Xpk - 1j*omeg**3*Xk**4*Xpk + 4*h*h[::-1,:]**2*rhok*Xpk**2 - 4.*1j*h*h[::-1,:]*omeg*R2*rhok*Xpk**2 - 4.*1j*h**2*omeg*rhok*Xk*Xpk**2 - 4.*1j*h*h[::-1,:]*omeg*rhok*Xk*Xpk**2 + 4.*1j*h[::-1,:]**2*omeg*rhok*Xk*Xpk**2 - 4*h*omeg**2*R2*rhok*Xk*Xpk**2 + 4*h[::-1,:]*omeg**2*R2*rhok*Xk*Xpk**2 + 8*h*omeg**2*rhok*Xk**2*Xpk**2 + 4*h[::-1,:]*omeg**2*rhok*Xk**2*Xpk**2 - 4.*1j*omeg**3*R2*rhok*Xk**2*Xpk**2 + 4.*1j*omeg**3*rhok*Xk**3*Xpk**2 + 4*domegR2*f*rhok*(h[::-1,:] - 1j*omeg*X)*((0,-1)*h + omeg*X)*(f - 1j*omeg*Xp) - 4*domegh[::-1,:]*rhok*(h + 1j*omeg*X)*(1j*f + omeg*Xp)*(f*R2 - f*Xk + h*Xpk + 1j*omeg*Xk*Xp) - 2*h*h[::-1,:]**2*R2*rhok*Xppk + 2*h*h[::-1,:]**2*rhok*Xk*Xppk + 4.*1j*h*h[::-1,:]*omeg*R2*rhok*Xk*Xppk - 2.*1j*h[::-1,:]**2*omeg*R2*rhok*Xk*Xppk - 4.*1j*h*h[::-1,:]*omeg*rhok*Xk**2*Xppk + 2.*1j*h[::-1,:]**2*omeg*rhok*Xk**2*Xppk + 2*h*omeg**2*R2*rhok*Xk**2*Xppk - 4*h[::-1,:]*omeg**2*R2*rhok*Xk**2*Xppk - 2*h*omeg**2*rhok*Xk**3*Xppk + 4*h[::-1,:]*omeg**2*rhok*Xk**3*Xppk + 2.*1j*omeg**3*R2*rhok*Xk**3*Xppk - 2.*1j*omeg**3*rhok*Xk**4*Xppk)
	den1 = -homeg[::-1,:]**4*homeg**2

	num2 = -(4*rhok*s12*(R2 - X)*(f**2*Xk - 2*f*homeg*Xpk + 1j*f*omeg*Xk*Xpk - 1j*homeg*omeg*Xpk**2 + domegh*f*(-1j*f + omeg*Xp)) - homeg**2*homeg[::-1,:]*s2*XPP)
	den2 = homeg[::-1,:]**2*homeg**3

	num3 = 4*f**2*rhok*Xk*(-(R2*s1) + homeg*s2 + s1*X) - 4*f*rhok*(2*homeg**2*s2 + s1*(homeg[::-1,:] + 1j*omeg*(R2 - X))*Xk + homeg*(-2*R2*s1 + 2*s1*Xk - 1j*omeg*s2*X))*Xpk - 4.*1j*rhok*(homeg*(1j*homeg[::-1,:]*s1 - omeg*R2*s1 + homeg*omeg*s2) + (homeg + homeg[::-1,:])*omeg*s1*X)*Xpk**2 + 4.*1j*domegh*rhok*(f + 1j*omeg*Xp)*(f*(R2*s1 - homeg*s2 - s1*X) + homeg[::-1,:]*s1*Xp) + homeg**2*s1*(-R2 + X)*XPP
	den3 = -homeg[::-1,:]*homeg**4
	

	XdynTemp=4.*qdim*(num1/den1+num2/den2+num3/den3)
	XdynNew= dot(dot(XdynTemp,wQ),wOmeg)*1./(2.*pi)
	
	Xerror=0.
	if abs(XdynNew.imag)>10**(-10.):
  		Xerror=XdynNew.imag
	return XdynNew.real,Xerror

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


