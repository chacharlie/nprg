from numpy import *

def weigthDicho(etaZlowT,etaZhighT,a):
	# Il semblerait que des valeurs de a proche de 0.4 (ou plus ?) soit ce qui marche le mieux : c'est-a-dire que ca ne vaut pas le coup d'essayer de faire mieux que la dichotomie naive
	# (en fait, c'est probablement aussi la facon dont l'implementation de cette nouvelle dichotomie est faite qui n'est pas ideale ??) 

	tLowT = max(size(where(gradient(abs(gradient(etaZlowT)))<0.)),5)
	tHighT = max(size(where(gradient(abs(gradient(etaZhighT)))<0.)),5)
	tsum = float(tLowT+tHighT)

	return a + (1.-2.*a)*tLowT/tsum, a + (1.-2.*a)*tHighT/tsum
