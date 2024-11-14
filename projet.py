#%% Imports

import numpy as np
import matplotlib.pyplot as plt

#%% Fonction aide
def readlist(file) : return list(map(float,file.readline().split()))
# Cette fonction permet simplement de lire une ligne d'un fichier et d'enregister chaque valeur
# séparé par un espace comme variable

def trace(x, y, xlabel, ylabel, titre, xlim=0, ylim=0, save=False, nom=None):
    fig, ax = plt.subplots(1, 1)

    ax.plot(x, y)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.suptitle(titre)

    if xlim != 0 and isinstance(xlim, list):
        ax.set_xlim(xlim)
    if ylim != 0 and isinstance(ylim, list):
        ax.set_ylim(ylim)

    plt.show()

    if save == True and nom != None:
        fig.savefig(nom, dpi=900)

#%% Création et conversion des données

VLAC = 750 #Volts
VSST = 790 #Volts
RSST = 33*1e-3 #mOhm
RHOLAC = 131e-6 # Ohm/m
RHORAIL = 18e-6 # Ohm/m

Times = []
X = []
V = [] #km/h-1 vitesse du train à déterminer avec le fichier marche_train.txt
Acc = [] #km/h-2 accélération du train à déterminer avec le fichier marche_train.txt
RLAC1 = [] # valeurs dépendante de x
RLAC2 = [] # valeurs dépendante de x
Rrail1 = [] # valeurs dépendante de x
Rrail2 = [] # valeurs dépendante de x
R1 = [] # Résistance équivalente pour la partie supérieure du schéma de Thévenin, traversée par le courant I1
R2 = [] # Résistance équivalente pour la partie inférieure du schéma de Thévenin, traversée par le courant I2
Req = [] # Résistance équivalente totale du schéma de Thévenin
PLAC = [] #dépend de x

alpha = 0 # angle de la pente du chemin
M = 70*1e3 #tonnes masse du train
A0 = 780 #N constante forces
A1 = 6.4*1e-3 #N/tonnes constante accélération
B0 = 0.0 # constante nulle ?
B1 = 0.14*3600/(1e3*1e3) #N/tonnes/(km•h-1) constante
C0 = 0.3634*(3600**2)/(1e3*1e6) #N/tonnes/(km•h-1)^2 constante inverse vitesse
C1 = 0.0

#%% Ajout des données du fichier

FICHIER = "marche_train.txt" # nom du fichier dans le dossier
file = open(FICHIER, 'r') # indique que nous ouvrons le fichier en lecture uniquement
readedlist = readlist(file)


while readedlist:
    Time, Position = readedlist
    Times.append(Time)
    X.append(Position)
    readedlist = readlist(file)

Times = np.array(Times)
X = np.array(X)

#%% Calcule vitesse, accélération, forces et puissances

V.append(0) #TODO revoir différence fini ordre supérieure
for i in range(len(X)-1):
    v = (X[i+1]-X[i])/(Times[i+1]-Times[i])
    V.append(v)

V = np.array(V)

Acc.append(0)
for i in range(len(X)-1):
    a = (V[i+1]-V[i])/(Times[i+1]-Times[i])
    Acc.append(a)

Acc = np.array(Acc)

FR = (A0 + A1*M) + (B0 + B1*M)*V + (C0 + C1*M)*V**2 # Force resistive

Fm = M*Acc + M*9.81*np.sin(alpha) + FR # Force mécanique - ici alpha = 0

Pm = Fm*V



#%% Partie électronique
#Calcul de RLAC1, RLAC2, Rrail1, Rrail2 en fonction de x
RLAC1.append(0)
RLAC2.append(0)
Rrail1.append(0)
Rrail2.append(0)

for i in range(len(X)-1) :
    rlac1 = RHOLAC*X[i]
    RLAC1.append(rlac1)

RLAC1 = np.array(RLAC1)

for i in range(len(X)-1) :
    rlac2 = RHOLAC*(2000 - X[i])
    RLAC2.append(rlac2)

RLAC2 = np.array(RLAC2)

for i in range(len(X)-1) :
    rrail1 = RHORAIL*X[i]
    Rrail1.append(rrail1)

Rrail1 = np.array(Rrail1)

for i in range(len(X)-1) :
    rrail2 = RHORAIL*(2000 - X[i])
    Rrail2.append(rrail2)

Rrail2 = np.array(Rrail2)

# Après simplification du schéma par le théorème de Thévenin calcul de R1, R2 et Req :
R1.append(0)
R2.append(0)
Req.append(0)

for i in range(len(X)-1) :
    r1 = RSST + RLAC1[i] + Rrail1[i]
    R1.append(r1)

R1 = np.array(R1)

for i in range(len(X)-1) :
    r2 = RSST + RLAC2[i] + Rrail2[i]
    R2.append(r2)

R2 = np.array(R2)

for i in range(len(X)-1) :
    req = (R1[i]*R2[i])/(R1[i]+R2[i])
    Req.append(req)

Req = np.array(Req)

# Calcule de PLAC :



Vtrain = (VSST + np.sqrt(VSST**2 - 4*Req*PLAC))/2

Itrain = (VSST - Vtrain)/REQ
I1 = Itrain * REQ/R1

PSST = VSST * I1

#%% Graphique

fig, ax = plt.subplots(3, 1)
ax[0].plot(Times, X)
ax[1].plot(Times, V)
ax[2].plot(Times, Pm*1e-6)
ax[2].plot(Times, np.zeros(len(Times)), '--', color='red')

ax[0].set_xlim([0,145])
ax[1].set_xlim([0,145])
ax[2].set_xlim([0,145])

ax[0].set_ylim([0,1300])
ax[2].set_ylim([-1,1.3])


ax[0].set_ylabel("x(t) [m]")
ax[1].set_ylabel("v(t) [m/s]")
ax[2].set_ylabel("Pm(t) [MW]")

ax[2].set_xlabel("t [s]")

fig.suptitle("Position, Vitesse et Puissance mécanique en fonction du temps")

plt.show()
fig.savefig("XVPm.pdf", dpi=900)

fig, ax = plt.subplots(3, 1)
ax[0].plot(Times, X)
ax[1].plot(Times, Pm*1e-6)
ax[1].plot(Times, np.zeros(len(Times)), '--', color='red')
ax[2].plot(Times, Vtrain)
ax[2].plot(Times, np.ones(len(Times))*500, '--', color='red')
ax[2].plot(Times, np.ones(len(Times))*790, '--', color='green')


ax[0].set_xlim([0,145])
ax[1].set_xlim([0,145])
ax[2].set_xlim([0,145])

ax[0].set_ylim([0,1300])
ax[1].set_ylim([-1,1.3])
ax[2].set_ylim([425,950])

ax[0].set_ylabel("x(t) [m]")
ax[1].set_ylabel("Pm(t) [MW]")
ax[2].set_ylabel("Vtrain(t) [V]")

ax[2].set_xlabel("t [s]")

fig.suptitle("Position, Puissance mécanique et Tension du train en fonction du temps")

plt.show()
fig.savefig("XpmVtrain.pdf", dpi=900)

trace(Times, Acc, "Temps [s]", "Accélération [$m/s^2$]", "Accélération en fonction du temps", [0, 140], [min(Acc[:140])-10,max(Acc[:140])+10], save=True, nom = "Acc.pdf")
