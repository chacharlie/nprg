from pylab import *
from numpy import *

# physique
dim=3
if dim==2 or dim==3 or dim==4:
	dims = str(dim)
else: 
	dims = str(int(dim*10))
NN=1

# numerique
NQ = 50         # nombre de pas pour les impulsions
Nomeg = 200     # nombre de pas pour les frequences
Nrho = 3*30        # nombre de pas pour le potentiel

# geometrique
LQ = 4.2                # taille du domaine selon q
Lrho = 1*0.0984*NN/(2**(-1-dim)*pi**((-dim/2))/math.gamma(dim/2)) # taille du domaine selon rho
Lrhos = "%.2f" % Lrho
drho = Lrho/Nrho        # pas de potentiel

# variables globales
rho = linspace(0,Lrho,Nrho)

#Runge-Kutta adaptatif
atol = 1e-6	# tolerance absolue sur l'erreur dans le Runge-Kutta
rtol = 0*.1e-6  	# tolerance relative

Ndicho=30
propDicho=0.5

beta=0
alpha=2.
kappa=3.5#1.25
#kappa=2.14251937166

RKadaptatif=True
NT=80000

model='A'
approx=2
choixRegu=1
choiceReguQ=0

edgeOrder5=True
diffOrder=5
edgeOrder=5
T=-30.

afficheNu=True
afficheEta=True
afficheVZX=True

printKappa=True

omegaRangeDico={'0.01': 1000,'0.02':500,'0.05':400,'0.07':142.857,'0.1':100.,'0.21':47.619,'0.25':40.,'0.5':20.,'0.75':13.3333,'1.0':10,'5.0':20,'10.0':10}
if choixRegu==1:
	Lomeg = 50.		# taille du domaine selon omega
elif choixRegu==2:
	Lomeg=omegaRangeDico[str(beta)]

#folderPath='results/N1d3alpha2NT40000Nrho30NQ50/'
#folderPath='results/N1d3alpha2Nrho'+str(Nrho)+'NQ'+str(NQ)+'/'
folderPath='results/March16/N'+str(NN)+'d'+dims+'Nrho'+str(Nrho)+'NQ'+str(NQ)+'/'

str1=model+'-'+str(approx)+'-'+str(Ndicho)+'-'
str3=str(Nomeg)+'-'+str(Lomeg)+'-'+str(beta)+'-'+str(alpha)+'-'+str(kappa)+'-'+str(choixRegu)+'-'+str(propDicho)+'-'+str(diffOrder)+'-'+str(edgeOrder)+'-Tmax'+str(T)+'-'+str(choiceReguQ)+'-'+Lrhos

if RKadaptatif:
	str2=str(atol)+'-'+str(rtol)+'-'
else:	
	str2=str(NT)+'-'

filePath=str1+str2+str3

fileName = folderPath + filePath+'testEta'# +'simpleLrho'#+'simpleInteg' #+ 'doubledLrho'# +'modifiedV0' #+'doubledLrho'# + 'test'

data=load(fileName+'.npz')
matrixEtaZ=data['etaZResults']
matrixEtaX=data['etaXResults']
matrixy=data['yResults']
#matrixy=data['Vresults']
matrixKappa=data['kappaResults']
matrixT=data['tResults']

