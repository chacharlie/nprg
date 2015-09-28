from pylab import *

from global_variables import model
if model=='ON':
	from eq_flow import *
elif model=='A':
	from A_eq_flow import *


def stepper(htry,y,dty,etaZ,etaX):
#Attempts a step with stepsize htry. On output, y and x are replaced by their new values, hdid
#is the stepsize that was actually accomplished, and hnext is the estimated next stepsize.

	
	h = htry # Set stepsize to the initial trial value
	control = ControlError(h)
	
	while 1==1:
		y_new,dty_new,etaZ_new,etaX_new,yerr = computeFlow(y,dty,etaZ,etaX,h) #Take a step.
		err = computeError(yerr,y,y_new) 	 # Evaluate accuracy of the step
		varSuccess = control.success(err)
		h = control.h
		if varSuccess:
			break
	
		#Step rejected. Try again with reduced h set by controller
		if (abs(h) <= 10**(-5)):
			print "stepsize underflow in StepperDopr5"
	
	dty = dty_new	 #Reuse last derivative evaluation for next step.
	y = y_new
	hnext = control.hnext

	return hnext,h,y_new,dty_new,etaZ_new,etaX_new 

def computeFlow(y,dty,etaZ,etaX,h):
	#Given values for Nrho variables V[0..Nrho-1] and their derivatives dtV[0..n-1], use the
	#fifth-order Dormand-Prince Runge-Kutta method to advance the solution over an interval h and
	#store the incremented variables in V_new[0..n-1]. Also store an estimate of the local truncation
	#error in Verr using the embedded fourth-order method.

	#constants for Runge-Kutta algorithm
	c2=0.2;c3=0.3;c4=0.8;c5=8.0/9.0;a21=0.2;a31=3.0/40.0;
	a32=9.0/40.0;a41=44.0/45.0;a42=-56.0/15.0;a43=32.0/9.0;a51=19372.0/6561.0;
	a52=-25360.0/2187.0;a53=64448.0/6561.0;a54=-212.0/729.0;a61=9017.0/3168.0;
	a62=-355.0/33.0;a63=46732.0/5247.0;a64=49.0/176.0;a65=-5103.0/18656.0;
	a71=35.0/384.0;a73=500.0/1113.0;a74=125.0/192.0;a75=-2187.0/6784.0;
	a76=11.0/84.0;e1=71.0/57600.0;e3=-71.0/16695.0;e4=71.0/1920.0;
	e5=-17253.0/339200.0;e6=22.0/525.0;e7=-1.0/40.0

	#First step.
	ytemp=y+h*a21*dty
	k2 = eqFlow(ytemp,etaZ,etaX,0)		# Second step.
	ytemp=y+h*(a31*dty+a32*k2)
	k3 = eqFlow(ytemp,etaZ,etaX,0)		# Third step.
	ytemp=y+h*(a41*dty+a42*k2+a43*k3)
	k4 = eqFlow(ytemp,etaZ,etaX,0)		# Fourth step.
	ytemp=y+h*(a51*dty+a52*k2+a53*k3+a54*k4)
	k5 = eqFlow(ytemp,etaZ,etaX,0)		# Fifth step.
	ytemp=y+h*(a61*dty+a62*k2+a63*k3+a64*k4+a65*k5)
	k6 = eqFlow(ytemp,etaZ,etaX,0) 		# Sixth step.

	# Now accumulate increments with proper weights.
	y_new = y+h*(a71*dty+a73*k3+a74*k4+a75*k5+a76*k6)
	dty_new,etaZ_new,etaX_new,yerror = eqFlow(y_new,etaZ,etaX,1) 	# Will also be first evaluation for next step.

	# Estimate error as difference between fourth- and fifth-order methods.
	yerr=h*(e1*dty+e3*k3+e4*k4+e5*k5+e6*k6+e7*dty_new)

	return y_new,dty_new,etaZ_new,etaX_new,yerr

def computeError(yerr,y,y_new):
	# Use Verr to compute norm of the scaled error estimate. A value less than one means the step was successful
	
	err=0.
	Nerr=size(y)
	for i in range(Nerr):		# a ameliorer : pas besoin de faire une boucle...
		sk = atol+rtol*max(abs(y[i]),abs(y_new[i]))
		err += (yerr[i]/sk)**2
	return (err/Nerr)**0.5

class ControlError:

	def __init__(self,hh = 0.):
		self.errold = 10**(-4.)
		self.reject = False
		self.hnext = 0.
		self.h = hh

	def success(self,err):
	#Returns true if err 1, false otherwise. If step was successful, sets hnext to the estimated
	#optimal stepsize for the next step. If the step failed, reduces h appropriately for another try.
	
		b=0.0;a=0.2-b*0.75;safe=0.9;minscale=0.2;maxscale=10.0;
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
			

