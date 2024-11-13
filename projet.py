#%% Imports
import numpy as np
import matplotlib.pyplot as plt
#%% Fonction aide
def readlist(file) : return list(map(float,file.readline().split()))
# Cette fonction permet simplement de lire une ligne d'un fichier et d'enregister chaque valeur
# séparé par un espace comme variable
#%% Création des données

VSST = 790 #Volts
RSST = 33*1e-3 #mOhm
RHOLAC = 131e-6 # Ohm/m
RHORAIL = 18e-6 # Ohm/m
PLAC = 35*1e3 #W

FICHIER = "marche_train.txt" # nom du fichier dans le dossier
file = open(FICHIER, 'r') # indique que nous ouvrons le fichier en lecture uniquement
readedlist = readlist(file)

Times = []
X = []

while readedlist:
    Time, Position = readedlist
    Times.append(Time)
    X.append(Position)
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
M = 70*1e3 #tonnes masse du train
A0 = 780 #N constante forces
A1 = 6.4*1e-3 #N/tonnes constante accélération
B0 = 0.0 # constante nulle ?
B1 = 0.14*3600/(1e3*1e3) #N/tonnes/(km•h-1) constante
C0 = 0.3634*(3600**2)/(1e3*1e6) #N/tonnes/(km•h-1)^2 constante inverse vitesse
C1 = 0.0

FR = (A0 + A1*M) + (B0 + B1*M)*V + (C0 + C1*M)*V**2 # Force resistive

Fm = M*Acc + M*9.81*np.sin(alpha) + FR # Force mécanique - ici alpha = 0

Pm = Fm*V

#%% Graphique
fig, ax = plt.subplots(3, 1)
ax[0].plot(Times, X)
ax[1].plot(Times, V)
ax[2].plot(Times, Pm*1e-6)
ax[2].plot(Times, np.zeros(len(Times)), '--', color='red')

ax[0].set_xlim([0,140])
ax[1].set_xlim([0,140])
ax[2].set_xlim([0,140])

ax[0].set_ylabel("x(t) [m]")
ax[1].set_ylabel("v(t) [m/s]")
ax[2].set_ylabel("Pm(t) [MW]")

ax[2].set_xlabel("t(s)")

fig.suptitle("Position, Vitesse et Puissance mécanique en fonction du temps")

plt.show()


#%% Partie électronique
R1 = RSST + (RHOLAC+RHORAIL)*2000 # TODO à vérifier
R1 = np.ones(len(Times))*R1
REQ = R1**2/(2*R1) # Car Somme de résistance en parallèle - 1/Req = 1/R1 + 1/R2
Vtrain = (VSST + np.sqrt(VSST**2 - 4*REQ*Pm))/2

Itrain = (VSST - Vtrain)/REQ
I1 = Itrain * REQ/R1
PSST = VSST * I1

fig, ax = plt.subplots(3, 1)
ax[0].plot(Times, X)
ax[1].plot(Times, Pm*1e-6)
ax[1].plot(Times, np.zeros(len(Times)), '--', color='red')
ax[2].plot(Times, Vtrain)
ax[2].plot(Times, np.ones(len(Times))*500, '--', color='red')
ax[2].plot(Times, np.ones(len(Times))*790, '--', color='green')


ax[0].set_xlim([0,140])
ax[1].set_xlim([0,140])
ax[2].set_xlim([0,140])

ax[0].set_ylim([0,1300])
ax[1].set_ylim([-0.8,1.3])
ax[2].set_ylim([425,950])

plt.show()
