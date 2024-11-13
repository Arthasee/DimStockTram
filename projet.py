#%% Imports
import numpy as np
import matplotlib.pyplot as plt
#%% Fonction aide
def readlist(file) : return list(map(float,file.readline().split()))
# Cette fonction permet simplement de lire une ligne d'un fichier et d'enregister chaque valeur
# séparé par un espace comme variable
#%% Création des données

VSST = 790 #Volts
RSST = 33 #mOhm
RHOLAC = 131e-6 # Ohm/m
RHORAIL = 18e-6 # Ohm/m

FICHIER = "marche_train.txt" # nom du fichier dans le dossier
file = open(FICHIER, 'r') # indique que nous ouvrons le fichier en lecture uniquement
readedlist = readlist(file)

Time, Position = readedlist

Times = []
X = []

while readedlist:
    Time, Position = readedlist
    Times.append(Time/3600.) # On convertit en heure
    X.append(Position*10**-3) # On convertit en km
    readedlist = readlist(file)

Times = np.array(Times)
X = np.array(X)

V = [] #km/h-1 vitesse du train à déterminer avec le fichier marche_train.txt
V.append(0)
for i in range(len(X)-1):
    v = (X[i+1]-X[i])/(Times[i+1]-Times[i])
    V.append(v)

V = np.array(V)

Acc = [] #km/h-2 accélération du train à déterminer avec le fichier marche_train.txt
Acc.append(0)
for i in range(len(X)-1):
    a = (V[i+1]-V[i])/(Times[i+1]-Times[i])
    Acc.append(a)

Acc = np.array(Acc)

alpha = 0 # angle de la pente du chemin
M = 70 #tonnes masse du train
A0 = 780 #N constante forces
A1 = 6.4 #N/tonnes constante accélération
B0 = 0.0 # constante nulle ?
B1 = 0.14 #N/tonnes/(km•h-1) constante
C0 = 0.3634 #N/tonnes/(km•h-1)^2 constante inverse vitesse
C1 = 0.0

FR = (A0 + A1*M) + (B0 + B1*M)*V + (C0 + C1*M)*V**2 # Force resistive

Fm = M*Acc + M*9.81*np.sin(alpha) + FR # Force mécanique - ici alpha = 0

Pm = Fm*V*10**-7

#%% Graphique
fig, ax = plt.subplots(3, 1)
ax[0].plot(Times*3600, X)
ax[1].plot(Times*3600, V)
ax[2].plot(Times*3600, Pm)

ax[0].set_xlim([0,140])
ax[1].set_xlim([0,140])
ax[2].set_xlim([0,140])

ax[0].set_ylim([0,1.5])

ax[2].set_xlabel("t(s)")

plt.show()

#%% Partie électronique
R1 = RSST*1e-3 + RHOLAC*X+RHORAIL*X # TODO à vérifier
REQ = R1/2 # Car Somme de résistance en parallèle - 1/Req = 1/R1 + 1/R2 - TODO à vérifier
Vtrain = 1./2. * (VSST + np.sqrt(VSST**2 - 4*REQ*(VSST**2/(4*REQ))))# TODO à modifier
