from pylab import *
from numpy import *

# physique
dim=3.                  # dimension d'espace
NN=1.                   # dimension des spins

# numerique
NT = 40000              # nombre de pas de temps (de RG)
NQ = 50         # nombre de pas pour les impulsions
Nomeg = 50     # nombre de pas pour les frequences
Nrho= 30                # nombre de pas pour le potentiel

# geometrique
T = -30.                                # taille du domaine selon t
LQ = 4.2                # taille du domaine selon q
Lomeg = 12.5              # taille du domaine selon omega
Lrho = 0.0984*NN/(2**(-1-dim)*pi**((-dim/2))/math.gamma(dim/2)) # taille du domaine selon rho
dt = T/NT                       # pas de temps
drho = Lrho/Nrho        # pas de potentiel

# variables globales
rho = linspace(0,Lrho,Nrho)


Ndicho=30
beta=0.1
kappa=4.5	


folderPath='results/N1d3alpha2NT40000Nrho30NQ50/'

fileName=folderPath+'Veta-'+str(Ndicho)+'-'+str(Nomeg)+'-'+str(Lomeg)+'-'+str(beta)+'-'+str(kappa)+'-moinsomega-flotXZcomplets'

data=load(fileName+'.npz')
matrixEtaZ=data['etaZResults']
matrixEtaX=data['etaXResults']
matrixV=data['Vresults']

indexPlotV=0
