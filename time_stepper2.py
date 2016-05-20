from pylab import *

from Dopr853_constants import *

from global_variables import model
if model=='ON':
	from eq_flow import *
elif model=='A':
	from A_eq_flow import *


def stepper(htry,y,dty):
#Attempts a step with stepsize htry. On output, y and x are replaced by their new values, hdid
#is the stepsize that was actually accomplished, and hnext is the estimated next stepsize.

	
	h = htry # Set stepsize to the initial trial value
	control = ControlError853(h)

	etaZ,etaX,rho0 = computeEta(y)

	while 1==1:
		y_new,yerr,yerr2 = computeFlow853(y,dty,etaZ,etaX,h) #Take a step.
	#	etaZ,etaX,rho0 = computeEta(y_new) #testEta
		err = computeError853(yerr,yerr2,y,y_new,h)
		varSuccess = control.success(err)
		h = control.h
		if varSuccess:
			break
	
		#Step rejected. Try again with reduced h set by controller
		if (abs(h) <= 10**(-5)):
			print "stepsize underflow in StepperDopr853, h=",h
	

	dty_new = eqFlow(y_new,etaZ,etaX,0)
	dty = dty_new
	y = y_new
	hnext = control.hnext

	return hnext,h,y_new,dty_new,etaZ,etaX,rho0 


def computeEta(y):
	# etaZ_new = A + B*etaZ + C*etaX
	_,A,Ap,rho0 = eqFlow(y,0.,0.,1)
	_,etaZ10,etaX10,_ = eqFlow(y,1.,0.,1)
	_,etaZ01,etaX01,_ = eqFlow(y,0.,1.,1)

	B,Bp = etaZ10-A, etaX10-Ap
	C,Cp = etaZ01-A, etaX01-Ap
	
	etaZ = (A*(1.-Cp)+C*A)/((1.-B)*(1.-Cp)-C*Bp)
	etaX = (Ap+Bp*etaZ)/(1.-Cp)

	return etaZ, etaX,rho0

def computeFlow853(y,dty,etaZ,etaX,h):
	ytemp=y+h*a21*dty
	k2 = eqFlow(ytemp,etaZ,etaX,0)
	
	ytemp=y+h*(a31*dty+a32*k2)
	k3 = eqFlow(ytemp,etaZ,etaX,0)
	
	ytemp=y+h*(a41*dty+a43*k3)
	k4 = eqFlow(ytemp,etaZ,etaX,0)
	
	ytemp=y+h*(a51*dty+a53*k3+a54*k4)
	k5 = eqFlow(ytemp,etaZ,etaX,0)
	
	ytemp=y+h*(a61*dty+a64*k4+a65*k5)
	k6 = eqFlow(ytemp,etaZ,etaX,0)
	
	ytemp=y+h*(a71*dty+a74*k4+a75*k5+a76*k6)
	k7 = eqFlow(ytemp,etaZ,etaX,0)
	
	ytemp=y+h*(a81*dty+a84*k4+a85*k5+a86*k6+a87*k7)
	k8 = eqFlow(ytemp,etaZ,etaX,0)
	
	ytemp=y+h*(a91*dty+a94*k4+a95*k5+a96*k6+a97*k7+a98*k8)
	k9 = eqFlow(ytemp,etaZ,etaX,0)
	
	ytemp=y+h*(a101*dty+a104*k4+a105*k5+a106*k6+a107*k7+a108*k8+a109*k9)
	k10 =eqFlow(ytemp,etaZ,etaX,0)
	
	ytemp=y+h*(a111*dty+a114*k4+a115*k5+a116*k6+a117*k7+a118*k8+a119*k9+a1110*k10)
	k2 = eqFlow(ytemp,etaZ,etaX,0)
	
	ytemp=y+h*(a121*dty+a124*k4+a125*k5+a126*k6+a127*k7+a128*k8+a129*k9+a1210*k10+a1211*k2)
	k3 = eqFlow(ytemp,etaZ,etaX,0)
	
	k4 = b1*dty+b6*k6+b7*k7+b8*k8+b9*k9+b10*k10+b11*k2+b12*k3
	y_new=y+h*k4

	yerr=k4-bhh1*dty-bhh2*k9-bhh3*k3
	yerr2=er1*dty+er6*k6+er7*k7+er8*k8+er9*k9+er10*k10+er11*k2+er12*k3

	return y_new,yerr,yerr2


def computeError853(yerr,yerr2,y,y_new,h):
	# Use Verr to compute norm of the scaled error estimate. A value less than one means the step was successful
	
	err=0.;err2=0.
	Nerr=size(y)
	for i in range(Nerr):		# a ameliorer : pas besoin de faire une boucle...
		sk = atol+rtol*max(abs(y[i]),abs(y_new[i]))
		err += (yerr[i]/sk)**2
		err2 += (yerr2[i]/sk)**2
	
	deno=err+0.01*err2
	if deno<=0.:
		deno=1.

	return abs(h)*err*(1./(Nerr*deno))**0.5

class ControlError853:

	def __init__(self,hh = 0.):
		self.errold = 10**(-4.)
		self.reject = False
		self.hnext = 0.
		self.h = hh

	def success(self,err):
	#Returns true if err 1, false otherwise. If step was successful, sets hnext to the estimated
	#optimal stepsize for the next step. If the step failed, reduces h appropriately for another try.
	
		b=0.0;a=1./8.-b*0.2;safe=0.9;minscale=0.333;maxscale=6.0;
		#Set beta to a nonzero value for PI control. beta=0.04 or 0.08 is a good default.
	
		if err<=1.0:	 # Step succeeded. Compute hnext.
			if err==0.0:
				scale=maxscale
			else:	# PI control if b != 0.
				scale=safe*err**(-a)*self.errold**b
				if scale<minscale:
					 scale=minscale # Ensure minscale <= hnext/h <= maxscale.
				if scale>maxscale:
					 scale=maxscale
			if self.reject:	# Dont let step increase if last one was rejected
				self.hnext=self.h*min(scale,1.0)
			else:
				self.hnext=self.h*scale
			self.errold = max(err,10**(-4.))	# Bookkeeping for next call.
			self.reject = False
			return True
	
		else:	# Truncation error too large, reduce stepsize.
			scale=max(safe*err**(-a),minscale)
			self.h *= scale
			self.reject = True
			return False
			

