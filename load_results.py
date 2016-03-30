from pylab import *
from numpy import *

# physique
dim=3
NN=1

# numerique
NQ = 50         # nombre de pas pour les impulsions
Nomeg = 200     # nombre de pas pour les frequences
Nrho= 30        # nombre de pas pour le potentiel

# geometrique
LQ = 4.2                # taille du domaine selon q
Lomeg = 50.              # taille du domaine selon omega
Lrho = 0.0984*NN/(2**(-1-dim)*pi**((-dim/2))/math.gamma(dim/2)) # taille du domaine selon rho
drho = Lrho/Nrho        # pas de potentiel

# variables globales
rho = linspace(0,Lrho,Nrho)

#Runge-Kutta adaptatif
atol = 1e-6	# tolerance absolue sur l'erreur dans le Runge-Kutta
rtol = 0*.1e-6  	# tolerance relative

Ndicho=30
propDicho=0.5

beta=0.2
alpha=2.
kappa=3.5
#kappa=2.14251937166

RKadaptatif=False
NT=40000

model='A'
approx=4
choixRegu=1

edgeOrder5=True
diffOrder=5
edgeOrder=5
T=-30.

afficheNu=True
afficheEta=True
afficheVZX=True

printKappa=True


#folderPath='results/N1d3alpha2NT40000Nrho30NQ50/'
#folderPath='results/N1d3alpha2Nrho'+str(Nrho)+'NQ'+str(NQ)+'/'
folderPath='results/March16/N'+str(NN)+'d'+str(dim)+'Nrho'+str(Nrho)+'NQ'+str(NQ)+'/'

str1=model+'-'+str(approx)+'-'+str(Ndicho)+'-'
str3=str(Nomeg)+'-'+str(Lomeg)+'-'+str(beta)+'-'+str(alpha)+'-'+str(kappa)+'-'+str(choixRegu)+'-'+str(propDicho)+'-'+str(diffOrder)+'-'+str(edgeOrder)+'-Tmax'+str(T)

if RKadaptatif:
	str2=str(atol)+'-'+str(rtol)+'-'
else:	
	str2=str(NT)+'-'

filePath=str1+str2+str3

#filePath=model+'-'+str(approx)+'-'+str(Ndicho)+'-'+str(atol)+'-'+str(rtol)+'-'+str(Nomeg)+'-'+str(Lomeg)\
#	+'-'+str(beta)+'-'+str(kappa)+'-'+str(choixRegu)+'-'+str(propDicho)+'-moinsomega-stepper2-diff-'+str(diffOrder)+'-'+str(edgeOrder)+'-Tmax'+str(T)+'testPlusDeBreak'

#filePath=model+'-'+str(approx)+'-'+str(Ndicho)+'-'+str(atol)+'-'+str(rtol)+'-'+str(Nomeg)+'-'+str(Lomeg)\
#	+'-'+str(beta)+'-'+str(kappa)+'-'+str(choixRegu)+'-'+str(propDicho)+'-moinsomega-stepper2-edge5-'+str(edgeOrder5)

#filePath=model+'-'+str(approx)+'-'+str(Ndicho)+'-'+str(atol)+'-'+str(rtol)+'-'+str(Nomeg)+'-'+str(Lomeg)\
#	+'-'+str(beta)+'-'+str(kappa)+'-'+str(choixRegu)+'-'+str(propDicho)+'-moinsomega'

#filePath=model+'-'+str(Ndicho)+'-'+str(atol)+'-'+str(rtol)+'-'+str(Nomeg)+'-'+str(Lomeg)\
#	+'-'+str(beta)+'-'+str(kappa)+'-'+str(choixRegu)+'-'+str(propDicho)+'-moinsomega'


#filePath='Veta-'+str(Ndicho)+'-'+str(atol)+'-'+str(rtol)+'-'+str(Nomeg)+'-'+str(Lomeg)+\
#	'-'+str(beta)+'-'+str(kappa)+'-'+str(choixRegu)+'-'+str(propDicho)+'-moinsomega'

#filePath='Veta-'+str(Ndicho)+'-'+str(Nomeg)+'-'+str(Lomeg)+'-'+str(beta)+'-'+str(kappa)+'-'+str(choixRegu)+'-moinsomega-flotXZcomplets'

fileName = folderPath + filePath

data=load(fileName+'.npz')
matrixEtaZ=data['etaZResults']
matrixEtaX=data['etaXResults']
matrixy=data['yResults']
#matrixy=data['Vresults']
matrixKappa=data['kappaResults']
matrixT=data['tResults']

