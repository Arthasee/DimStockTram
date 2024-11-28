#%% Imports

import numpy as np
import matplotlib.pyplot as plt
from help_function import *

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
RLAC1 = [] # Résistance de la LAC entre la sous-station 1 et le train (valeurs dépendante de x)
RLAC2 = [] # Résistance de la LAC entre la sous-station 2 et le train (valeurs dépendante de x)
Rrail1 = [] # Résistance du rail entre la sous-station 1 et le train (valeurs dépendante de x)
Rrail2 = [] # Résistance du rail entre la sous-station 2 et le train (valeurs dépendante de x)
R1 = [] # Résistance équivalente pour la partie supérieure du schéma de Thévenin, traversée par le courant I1 (dépend de x)
R2 = [] # Résistance équivalente pour la partie inférieure du schéma de Thévenin, traversée par le courant I2 (dépend de x)
Req = [] # Résistance équivalente totale du schéma de Thévenin (dépend de x)
PLAC = [] #Puissance de la LAC (dépend de x)
Vtrain = [] #Tension du train à tout instant (dépend de x)
Itrain = [] #Intensité aux bornes du train à tout moment (dépend de x)
I1 = [] #Intensité de la partie supérieure du schéma de Thévenin
I2 = [] #Intensité de la partie inférieure du schéma de Thévenin

alpha = 0 # angle de la pente du chemin
M = 70*1e3 #tonnes masse du train
A0 = 780 #N constante forces
A1 = 6.4*1e-3 #N/tonnes constante accélération
B0 = 0.0 # constante nulle ?
B1 = 0.14*3600/(1e3*1e3) #N/tonnes/(km•h-1) constante
C0 = 0.3634*(3600**2)/(1e6) #N/(km•h-1)^2 constante inverse vitesse
C1 = 0.0

#%% Ajout des données du fichier

FICHIER = "marche_train.txt" # nom du fichier dans le dossier

Times, X = get_T_X(FICHIER)

#%% Calcule vitesse, accélération, forces et puissances

V = get_V(Times, X)

Acc = get_Acc(Times, V)

FR = (A0 + A1*M) + (B0 + B1*M)*V + (C0 + C1*M)*V**2 # Force resistive

Fm = M*Acc + M*9.81*np.sin(alpha) + FR # Force mécanique - ici alpha = 0

Pm = Fm*V

#%% Partie électronique
#Calcul de RLAC1, RLAC2, Rrail1, Rrail2 en fonction de x

RLAC1 = RHOLAC*X

RLAC2 = (X[-1] - X)*RHOLAC

Rrail1 = RHORAIL*X

Rrail2 = (X[-1] - X)*RHORAIL

# Après simplification du schéma par le théorème de Thévenin calcul de R1, R2 et Req :
R1 = RSST + RLAC1 + Rrail1

R2 = RSST + RLAC2 + Rrail2

Req = (R1*R2)/(R1+R2)

# Calcul de PLAC :
PLAC = VLAC**2/(RLAC1+RLAC2)

# Calcul de Vtrain :
for i in range(len(X)):
    if Pm[i]<0:
        racine = VSST**2 - 4*Req[i]*(Pm[i]*0.8)
    else:
        racine = VSST**2 - 4*Req[i]*(Pm[i]/0.8)
    racine = max(racine, 0)
    vtrain = (VSST + np.sqrt(racine))/2
    Vtrain.append(vtrain)

Vtrain = np.array(Vtrain)
#ici un test
# Calcul de Itrain :
Itrain = VSST - Vtrain/Req

# Calcul de I1 : 
# On sait d'après la loi des mailles que V1 - V2 = 0, donc V1 = V2, donc R1 * I1 = R2 * I2, donc I2 = (R1 * I1)/R2
# D'après la loi des noeuds, I1 + I2 = Itrain, donc en remplaçant I2 par son expression en fonction de I1 on obtient :
# I1 + (R1 * I1)/R2 = Itrain, donc I1(R2 + R1)/R2 = Itrain, donc I1 = (R2 * Itrain)/(R1 + R2)
I1 = (R2*Itrain)/(R1+R2)

# Calcul de I2 :
I2 = (R1*I1)/R2

# Calcul de la puissance de chaque sous-station : Psst = Vsst*Isst = Vsst**2 / Rsst
PSST = VSST**2 / RSST

#%% Graphique

# fig, ax = plt.subplots(3, 1)
# ax[0].plot(Times, X)
# ax[1].plot(Times, V)
# ax[2].plot(Times, Pm*1e-6)
# ax[2].plot(Times, np.zeros(len(Times)), '--', color='red')

# ax[0].set_xlim([0,145])
# ax[1].set_xlim([0,145])
# ax[2].set_xlim([0,145])

# ax[0].set_ylim([0,1300])
# ax[2].set_ylim([-1,1.3])


# ax[0].set_ylabel("x(t) [m]")
# ax[1].set_ylabel("v(t) [m/s]")
# ax[2].set_ylabel("Pm(t) [MW]")

# ax[2].set_xlabel("t [s]")

# fig.suptitle("Position, Vitesse et Puissance mécanique en fonction du temps")

# plt.show()
# fig.savefig("XVPm.pdf", dpi=900)

# fig, ax = plt.subplots(3, 1)
# ax[0].plot(Times, X)
# ax[1].plot(Times, Pm*1e-6)
# ax[1].plot(Times, np.zeros(len(Times)), '--', color='red')
# ax[2].plot(Times, Vtrain)
# ax[2].plot(Times, np.ones(len(Times))*500, '--', color='red')
# ax[2].plot(Times, np.ones(len(Times))*790, '--', color='green')


# ax[0].set_xlim([0,145])
# ax[1].set_xlim([0,145])
# ax[2].set_xlim([0,145])

# ax[0].set_ylim([0,1300])
# ax[1].set_ylim([-1,1.3])
# ax[2].set_ylim([425,950])

# ax[0].set_ylabel("x(t) [m]")
# ax[1].set_ylabel("Pm(t) [MW]")
# ax[2].set_ylabel("Vtrain(t) [V]")

# ax[2].set_xlabel("t [s]")

# fig.suptitle("Position, Puissance mécanique et Tension du train en fonction du temps")

# plt.show()
# fig.savefig("XpmVtrain.pdf", dpi=900)

# fig, ax = plt.subplots(2, 1)
# ax[0].plot(Times, Pm)
# # ax[1].plot(Times, PLAC)

# fig.suptitle("Les Puissances en fonction du temps")
# plt.show()

# fig, ax = plt.subplots(3, 1)
# ax[0].plot(Times, R1)
# ax[1].plot(Times, R2)
# ax[2].plot(Times, Req)

# fig.suptitle("Les Résistances R1, R2 et Req en fonction du temps")
# plt.show()

# fig, ax = plt.subplots(3, 1)
# ax[0].plot(Times, I1)
# ax[1].plot(Times, I2)
# ax[2].plot(Times, Itrain)

# fig.suptitle("Les Intensités I1, I2 et Itrain en fonction du temps")
# plt.show()

# trace(Times, Acc, "Temps [s]", "Accélération [$m/s^2$]", "Accélération en fonction du temps", [0, 140], [-2,2], save=True, nom = "Acc.pdf")
# trace(Times, PLAC, "Temps [s]", "PLAC", "PLAC en fonction du temps", save=True, nom = "PLAC.pdf")

#%% Ajout de la batterie (en sst inversible)

# Initialisation
Pbatt = np.zeros(len(Pm))
Ebatt = np.zeros(len(Pm))
Ebatt0 = 172*1e3*3/4
EbattMAX = 172*1e3
Ebatt[0] = Ebatt0
Prheos = np.zeros(len(Pm))
VtrainBatt = np.zeros(len(Pm))
Pelec = np.zeros(len(Pm))

# Remplissage de Pelec, rendement de perte 0.8 
for i in range(1, len(Pm)):
    if Pm[i]<=0:
        Pelec[i] = Pm[i]*0.8
    else:
        Pelec[i] = Pm[i]/0.8

    seuil = 0.5 * np.max(Pelec)

    # Loi de gestion de la batterie
    if Pelec[i]<0 and Ebatt[i-1] < EbattMAX:
        Pbatt[i] = abs(Pelec[i])
        Ebatt[i] = Ebatt[i-1] + Pbatt[i]*1/3600
        if Ebatt[i] > EbattMAX:
            Ebatt[i] = EbattMAX
            Prheos[i] = Prheos[i-1] + Pelec[i] + (Ebatt[i]-EbattMAX)*3600
    elif Pelec[i]> seuil and Pbatt[i-1] > 0:
        Pbatt[i] = Ebatt[i-1] * 3600
        if Pelec[i] > Pbatt[i]:
            Pelec[i] -= Pbatt[i]
        else:
            Ebatt[i] = Ebatt[i-1] - Pelec[i]*1/3600
            Pelec[i] = 0
    else:
        Ebatt[i] = Ebatt[i-1]
        Pbatt[i] = Pbatt[i-1]

    PLAC[i] = Pelec[i] - Pbatt[i] + Prheos[i]
    if PLAC[i] < 0:
        PLAC[i] = 0
    print("Pbatt =", Pbatt[i], "Ebatt =", Ebatt[i], "PLAC =", PLAC[i])
    racine = VSST**2 - 4*Req[i]*(PLAC[i]/0.8)
    racine = max(racine,0)
    vtrain = (VSST + np.sqrt(racine))/2
    VtrainBatt[i] = vtrain

# Affichage des solutions 
trace(Times, Ebatt, "Temps[s]", "Energie de la batterie", "Energie de la batterie en fonction du temps")
trace(Times, PLAC, "Temps[s]", "PLAC", "PLAC avec batterie en fonction du temps")
trace(Times, Pbatt, "Temps[s]", "puissance batterie", "puissance batterie en fonction du temps")
trace(Times, VtrainBatt, "Temps[s]", "Vtrain", "Vtrain avec batterie en fonction du temps") #, [0, 140]


#%% Dimmensionnement du système de stockage
#construction de l’ensemble des solutions non dominées pour les critères « Capacité en énergie de la batterie » et « Chute de tension maximale » (qui est la différence entre la tension nominale (750V) et la tension réelle mesurée aux bornes du train.)

# Méthode de Monte-Carlo

# Nombre de simulations Monte-Carlo
nbre_simulations = 1000 # fixé à priori

# Stockage des résultats
capacite_batterie = np.random.uniform(50, 200, nbre_simulations)  # Capacité de la batterie (en kWh)
chute_tension = np.random.uniform(10, 250, nbre_simulations)  # Chute de tension maximale (en V)

# Construire les solutions non dominées
solutions_non_dominees = []
for i in range(nbre_simulations):
    is_dominated = False
    for j in range(nbre_simulations):
        if (capacite_batterie[j] <= capacite_batterie[i] and chute_tension[j] <= chute_tension[i] and 
            (capacite_batterie[j] < capacite_batterie[i] or chute_tension[j] < chute_tension[i])):
            is_dominated = True
            break
    if not is_dominated:
        solutions_non_dominees.append(i)

# Affichage des solutions 
plt.scatter(capacite_batterie, chute_tension, color = 'skyblue', label='Ensemble des solutions par la méthode de Monté - Carlo')
plt.scatter(capacite_batterie[solutions_non_dominees], chute_tension[solutions_non_dominees], color='red', label='Solutions non dominées')
plt.xlabel('Capacité en énergie de la batterie (kWh)')
plt.ylabel('Chute de tension maximale (V)')
plt.legend()
plt.show()


# Méthode NGSA2 ( non-sorted algorithm system)

""" Notes: Objectif: trouver un enssemble de solutions non-dominées appelées Paréto optimales. Identification de front de paréto 
-> ensemble de solutions où aucune solution ne peut être améliorée dans un objectif sans détériorer un autre objectif.
1) Population N choisis aléatoirement
2) Individu --> gène , enssemble des gène de l'inidividu --> chromosome
3) À chaque génération (nombre de génération choisi préalablement) : 
- Tri des individus, premier front contient les solutions non dominées( Paréto optimales), le second front contient les solutions dominées uniquement par celles du premier front, et ainsi de suite.
50 pourcent des meilleurs individus seront choisis (dans chaque front de paréto j'imagine)
- Mesure de diversité, distance de regroupement
- Opération génétique classique vu en IA TP4: Mutations et Croisement pour éxplorer d'autres espaces de recherches
4) Àrret à la fin du cycle de génération
"""


