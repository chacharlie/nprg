from pylab import *
from numpy import *

# physique
dim=3.                  # dimension d'espace
NN=1.                   # dimension des spins

# numerique
NQ = 50         # nombre de pas pour les impulsions
Nomeg = 200     # nombre de pas pour les frequences
Nrho= 30                # nombre de pas pour le potentiel

# geometrique
LQ = 4.2                # taille du domaine selon q
Lomeg = 50              # taille du domaine selon omega
Lrho = 0.0984*NN/(2**(-1-dim)*pi**((-dim/2))/math.gamma(dim/2)) # taille du domaine selon rho
drho = Lrho/Nrho        # pas de potentiel

# variables globales
rho = linspace(0,Lrho,Nrho)

#Runge-Kutta adaptatif
atol = 10**(-6.)	# tolerance absolue sur l'erreur dans le Runge-Kutta
rtol = 0.*10**(-6.)  	# tolerance relative

Ndicho=55
propDicho=0.4

beta=0.
kappa=1.5	

model='A'
approx=4
choixRegu=1

afficheNu=True
indexPlotV=54


#folderPath='results/N1d3alpha2NT40000Nrho30NQ50/'
folderPath='results/N1d3alpha2Nrho30NQ50/'

fileName=folderPath+model+'-'+str(approx)+'-'+str(Ndicho)+'-'+str(atol)+'-'+str(rtol)+'-'+str(Nomeg)+'-'+str(Lomeg)\
	+'-'+str(beta)+'-'+str(kappa)+'-'+str(choixRegu)+'-'+str(propDicho)+'-moinsomega'

#fileName=folderPath+model+'-'+str(Ndicho)+'-'+str(atol)+'-'+str(rtol)+'-'+str(Nomeg)+'-'+str(Lomeg)\
#	+'-'+str(beta)+'-'+str(kappa)+'-'+str(choixRegu)+'-'+str(propDicho)+'-moinsomega'


#fileName=folderPath+'Veta-'+str(Ndicho)+'-'+str(atol)+'-'+str(rtol)+'-'+str(Nomeg)+'-'+str(Lomeg)+\
#	'-'+str(beta)+'-'+str(kappa)+'-'+str(choixRegu)+'-'+str(propDicho)+'-moinsomega'

#fileName=folderPath+'Veta-'+str(Ndicho)+'-'+str(Nomeg)+'-'+str(Lomeg)+'-'+str(beta)+'-'+str(kappa)+'-'+str(choixRegu)+'-moinsomega-flotXZcomplets'

data=load(fileName+'.npz')
matrixEtaZ=data['etaZResults']
matrixEtaX=data['etaXResults']
matrixy=data['yResults']
#matrixy=data['Vresults']
matrixKappa=data['kappaResults']
matrixT=data['tResults']



