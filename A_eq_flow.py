from pylab import *
from scipy.interpolate import interp1d
import numpy as np

from global_variables import *
from diff_op import *
from regu import *

def eqFlow(VZX,etaZ,etaX,computeEta):
	V=VZX[0:Nrho]
	Vp=d_rho(V)
	Vpp=d2_rho(V)
	if approx==1:
		X = 1.*ones((rho.size))
#		X = 1.*ones((rho.size))+alpha*beta # modifX
#		X = DD*ones((rho.size)) # modifX2
		Xp = zeros((rho.size))
		Xpp = zeros((rho.size))
		Z = 1.*ones((rho.size))
		Zp = zeros((rho.size))
		Zpp = zeros((rho.size))
	if approx==2:
		X = 1.*ones((rho.size))
		Xp = zeros((rho.size))
		Xpp = zeros((rho.size))
		Z=VZX[Nrho:2*Nrho]
		Zp=d_rho(Z)
		Zpp=d2_rho(Z)
	if approx==3:
		Z = 1.*ones((rho.size))
		Zp = zeros((rho.size))
		Zpp = zeros((rho.size))
		X=VZX[2*Nrho:3*Nrho]
		Xp=d_rho(X)
		Xpp=d2_rho(X)
	if approx==4:
		Z=VZX[Nrho:2*Nrho]
		Zp=d_rho(Z)
		Zpp=d2_rho(Z)
		X=VZX[2*Nrho:3*Nrho]
		Xp=d_rho(X)
		Xpp=d2_rho(X)

	Vdyn = zeros((rho.size))
	Verror= zeros((rho.size))
	Xdyn = zeros((rho.size))
	Xerror= zeros((rho.size))
	Zdyn = zeros((rho.size))
	Zerror= zeros((rho.size))


	if beta==0 and exact==True:
		s = -(etaZ*regu+2.*qq*regup)
	else:
		s1 = -(etaZ*qR1+qdqR1+(2.-etaZ+etaX)*qomegdomegR1)
#		s12 = -(etaZ*qR1[::-1,:]+qdqR1[::-1,:]+(2.-etaZ+etaX)*qomegdomegR12)#s1[::-1,:] 
		s12 = -(etaZ*qR1[::-1,:]+qdqR1[::-1,:]-(2.-etaZ+etaX)*qomegdomegR12) # r2s12modif
		s2 = -(etaX*R2+qdqR2+(2-etaZ+etaX)*omegdomegR2)
#		s2 = -((etaX+2.*aa)*R2+qdqR2+(2-etaZ+etaX)*omegdomegR2)	# r2modif
	
	for k in range(Nrho):
		if beta==0 and exact==True:
			h = qq*(regu+Z[k])+V[k]+2.*rho[k]*Vp[k]
			f= qq*Zp[k]+3.*Vp[k]+2.*rho[k]*Vpp[k]
			Vdyn[k],Verror[k] = eqVbeta0(f,s,h)
			Zdyn[k],Zerror[k] = eqZbeta0(f,s,h,Z[k],Zp[k],Zpp[k],rho[k])
			Xdyn[k],Xerror[k] = eqXbeta0(f,s,h,X[k],Xp[k],Xpp[k],rho[k])
		else:	
			h = qq*(R1+Z[k])+V[k]+2.*rho[k]*Vp[k]
			homeg = h+1j*omeg*X[k]
			f= qq*Zp[k]+3.*Vp[k]+2.*rho[k]*Vpp[k]
			g = homeg*Xp[k]-X[k]*f
			R2X=R2-X[k]
			Vdyn[k],Verror[k] = eqV(f,g,s1,s12,s2,homeg)
			Zdyn[k],Zerror[k] = eqZ(rho[k],X[k],Xp[k],Z[k],Zp[k],Zpp[k],R2X,s1,s12,s2,f,homeg)
			Xdyn[k],Xerror[k] = eqX(rho[k],X[k],Xp[k],Xpp[k],R2X,s1,s12,s2,f,h,homeg)
		
 	Vdim = (-2.+etaZ)*V+(-2.+dim+etaZ)*rho*Vp 
	Zdim = etaZ*Z+(-2.+dim+etaZ)*rho*Zp
	Xdim = etaX*X+(-2.+dim+etaZ)*rho*Xp

	dtV = Vdim+Vdyn
	dtZ = Zdim+Zdyn
	dtX = Xdim+Xdyn
	
	dtVZX = np.array([dtV,dtZ,dtX])

#	VZXerror = [max(abs(Verror))]+[max(abs(Xerror))]+[max(abs(Zerror))]

	if computeEta==1:	
		dtZ0,dtX0,rho0 = findZX0(V,Vp,Zp,Xp,dtV,dtZ,dtX)
		return dtVZX.flatten(),etaZ-dtZ0,etaX-dtX0,rho0 #,np.array(VZXerror)
	else:
		return dtVZX.flatten()

def eqVbeta0(f,s,h):
# le 2*q provient du chgmt de variable qd on integre sur q au lieu de y
	VdynTemp = -2.*q*(q**dim)*f*s/h**2
	VdynNew= dot(VdynTemp,wQ)
	return VdynNew,0

def eqZbeta0(f,s,h,z,zp,zpp,rhok):
	hp = z + regu + qq*regup	
	dqqh = 2.*regup+qq*regupp

	t1 = 2*rhok*f**2/h**2*(-hp+4./dim*qq*hp**2/h-2./dim*qq*dqqh)
	t2 = 4.*rhok*zp*f/h*(1.-2./dim*qq*hp/h) 
	t3 = 2./dim*rhok*zp**2*qq/h
	t4 = -0.5*zp-rhok*zpp

	ZdynTemp = 2.*q*2.*(q**dim)*s/h**2*(t1+t2+t3+t4)
#	ZdynTemp = 2.*q*(q**dim)*s/h**2*(t1+t2+t3+t4) # factor2
	ZdynNew= dot(ZdynTemp,wQ)
	return ZdynNew,0

def eqXbeta0(f,s,h,x,xp,xpp,rhok):
	t1 = 3./2.*f**2/h**2*x*rhok
	t2 = -4.*rhok*xp*f/h
	t3 = 0.5*xp+rhok*xpp

	XdynTemp = -2.*q*2.*(q**dim)*s/h**2*(t1+t2+t3) 
#	XdynTemp = -2.*q*(q**dim)*s/h**2*(t1+t2+t3) #factor2
	XdynNew= dot(XdynTemp,wQ)
	return XdynNew,0


def eqV(f,g,s1,s12,s2,homeg):
	term1 = (g[::-1,:] + f*R2)*s1/(homeg[::-1,:]*homeg**2)
	term2 = (g + f*R2)*s12/(homeg[::-1,:]**2*homeg)
	term3 = -f*s2/(homeg[::-1,:]*homeg)
  
  
	VdynTemp=4.*qdim*(term1+term2+term3)
	VdynNew= dot(dot(VdynTemp,wQ),wOmeg)*1./(2.*pi)
  
	Verror=0.
	if abs(VdynNew.imag)>1e-8:
		Verror=VdynNew.imag
		print 'Imaginary part of V: {0}'.format(Verror)
	return VdynNew.real,Verror


def eqZ(rhok,Xk,Xpk,Zk,Zpk,Zppk,R2X,s1,s12,s2,f,homeg):
	dqh = Zk + dqhnoZ	
	ZPP = Zpk + 2*rhok*Zppk
	fXp = f + 1j*omeg*Xpk
	s22 = s1*Xk+s2*homeg-R2*s1
	
	nume1 = 4*rhok*s12*R2X*(-(dim*homeg*(-(dqh*f*fXp) + (f + fXp)*homeg*Zpk)) + qq*(-4*dqh**2*f*fXp + 2*dqh*(f + fXp)*homeg*Zpk + homeg*(2*dqqh*f*fXp - homeg*Zpk**2))) - dim*homeg**3*homeg[::-1,:]*s2*ZPP
	deno1 = homeg**4*homeg[::-1,:]**2

	nume2 = 4*homeg[::-1,:]*rhok*s1*Xpk*(dim*homeg*(dqh*fXp - homeg*Zpk) + 2*qq*(-2*dqh**2*fXp + dqqh*fXp*homeg + dqh*homeg*Zpk)) - 4*rhok*s22*qq*(-4*dqh**2*f*fXp + 2*dqh*(f + fXp)*homeg*Zpk + homeg*(2*dqqh*f*fXp - homeg*Zpk**2)) - dim*homeg*(4*dqh*f*fXp*rhok*s22 - homeg*(4*f*rhok*s22*Zpk + 4*fXp*rhok*s22*Zpk + homeg*s1*R2X*ZPP))
	deno2 = homeg**5*homeg[::-1,:]
	
	nume3 = -(s12*(4*rhok*qq*(4*dqh[::-1,:]**2*fXp[::-1,:]*homeg**2*(f*R2X + homeg*Xpk) - 2*dqh[::-1,:]*homeg*homeg[::-1,:]*(2*dqR2*f*fXp[::-1,:]*homeg - 2*dqh*f*fXp[::-1,:]*R2X + homeg*((f + fXp[::-1,:])*R2X + homeg*Xpk)*Zpk) + homeg[::-1,:]*(4*dqh**2*f*homeg[::-1,:]*(fXp[::-1,:]*R2X + homeg[::-1,:]*Xpk) - 2*dqh*homeg*homeg[::-1,:]*(2*dqR2*f*fXp[::-1,:] + ((f + fXp[::-1,:])*R2X + homeg[::-1,:]*Xpk)*Zpk) + homeg*(2*dqqR2*f*fXp[::-1,:]*homeg*homeg[::-1,:] - 2*dqqh[::-1,:]*fXp[::-1,:]*homeg*(f*R2X + homeg*Xpk) + homeg[::-1,:]*(-2*dqqh*f*(fXp[::-1,:]*R2X + homeg[::-1,:]*Xpk) + homeg*Zpk*(2*dqR2*(f + fXp[::-1,:]) + R2X*Zpk))))) + dim*homeg*homeg[::-1,:]*(4*dqR2*f*fXp[::-1,:]*homeg*homeg[::-1,:]*rhok - 4*dqh[::-1,:]*fXp[::-1,:]*homeg*rhok*(f*R2X + homeg*Xpk) + homeg[::-1,:]*(-4*dqh*f*rhok*(fXp[::-1,:]*R2X + homeg[::-1,:]*Xpk) + 4*homeg*rhok*((f + fXp[::-1,:])*R2X + (homeg + homeg[::-1,:])*Xpk)*Zpk - homeg*homeg[::-1,:]*R2X*ZPP))))
	deno3 = homeg[::-1,:]**5*homeg**3
	
	ZdynTemp=4.*qdim*1./dim*(nume1/deno1+nume2/deno2+nume3/deno3)
	ZdynNew= dot(dot(ZdynTemp,wQ),wOmeg)*1./(2.*pi)
	
	Zerror=0.
	if abs(ZdynNew.imag)>1e-8:
		Zerror=ZdynNew.imag
		print 'Imaginary part of Z: {0}'.format(Zerror)
	return ZdynNew.real,Zerror
	
def eqX(rhok,Xk,Xpk,Xppk,R2X,s1,s12,s2,f,h,homeg):
	XPP = Xpk + 2*rhok*Xppk
	domegh=qdomegR1

	num1 = 4*homeg*rhok*s12*(-1j*domegh[::-1,:] + Xk)*(f*R2X + homeg*Xpk)*(f - 1j*omeg*Xpk) - homeg[::-1,:]*s12*(4*rhok*(1j*f + omeg*Xpk)*(domegR2*f*homeg - f*R2X*(domegh + 1j*Xk) + 1j*homeg*R2X*Xpk) + homeg[::-1,:]*(4*f*rhok*(-1j*domegh + Xk)*Xpk + homeg*Xpk*(R2X - 4*rhok*Xpk) + 2*homeg*R2X*rhok*Xppk))
	den1 = -homeg[::-1,:]**4*homeg**2

	num2 = (4*rhok*s12*R2X*(f**2*Xk - 2*f*homeg*Xpk + 1j*f*omeg*Xk*Xpk - 1j*homeg*omeg*Xpk**2 + domegh*f*(-1j*f + omeg*Xpk)) - homeg**2*homeg[::-1,:]*s2*XPP)
	den2 = homeg[::-1,:]**2*homeg**3

	num3 = 4*f**2*rhok*Xk*(-(R2*s1) + homeg*s2 + s1*Xk) - 4*f*rhok*(2*homeg**2*s2 + s1*(homeg[::-1,:] + 1j*omeg*R2X)*Xk + homeg*(-2*R2*s1 + 2*s1*Xk - 1j*omeg*s2*Xk))*Xpk - 4.*1j*rhok*(homeg*(1j*homeg[::-1,:]*s1 - omeg*R2*s1 + homeg*omeg*s2) + (homeg + homeg[::-1,:])*omeg*s1*Xk)*Xpk**2 + 4.*1j*domegh*rhok*(f + 1j*omeg*Xpk)*(f*(R2*s1 - homeg*s2 - s1*Xk) + homeg[::-1,:]*s1*Xpk) - homeg**2*s1*R2X*XPP
	den3 = -homeg[::-1,:]*homeg**4

	XdynTemp=4.*qdim*(num1/den1+num2/den2+num3/den3)
	XdynNew= dot(dot(XdynTemp,wQ),wOmeg)*1./(2.*pi)

	Xerror=0.
	if abs(XdynNew.imag)>1e-8:
  		Xerror=XdynNew.imag
		print 'Imaginary part of X: {0}'.format(Xerror)
	return XdynNew.real,Xerror

def findZX0(V,Vp,Zp,Xp,dtV,dtZ,dtX):
	if int(abs(sum(sign(V))))!=Nrho:
		VinterpRho0=interp1d(V.real,rho)
		rho0=float(VinterpRho0(0.))
		dtZinterp=interp1d(rho,dtZ)
		dtXinterp=interp1d(rho,dtX)
		dtVinterp=interp1d(rho,dtV)
		VpInterp=interp1d(rho,Vp)
		ZpInterp=interp1d(rho,Zp)
		XpInterp=interp1d(rho,Xp)
		dtRho0 = - dtVinterp(rho0)/VpInterp(rho0)
		return dtZinterp(rho0)+ZpInterp(rho0)*dtRho0,dtXinterp(rho0)+XpInterp(rho0)*dtRho0,rho0
	else:
		return dtZ[0]-Zp[0]*dtV[0]/Vp[0],dtX[0]-Xp[0]*dtV[0]/Vp[0],0.
