#%% Imports

import numpy as np
import matplotlib.pyplot as plt
from help_function import *
import copy
import random

#%% Création et conversion des données
NB_SIMU= 100
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
for ij in range(len(X)):
    if Pm[ij]<0:
        rac = VSST**2 - 4*Req[ij]*(Pm[ij]*0.8)
    else:
        rac = VSST**2 - 4*Req[ij]*(Pm[ij]/0.8)
    rac = max(rac, 0)
    vt = (VSST + np.sqrt(rac))/2
    Vtrain.append(vt)

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

def remplissage_P_elec(Pmi):
    # Initialisation
    pelec = np.zeros(len(Pmi))

    # Remplissage de Pelec, rendement de perte 0.8
    for i in range(1, len(Pmi)):
        if Pmi[i]<=0:
            pelec[i] = Pmi[i]*0.8
        else:
            pelec[i] = Pmi[i]/0.8
    print(np.max(pelec))
    return pelec

def gestion_batterie(pelec,EbattMAX,seuil):
    # Initialisation
    Pbatt = np.zeros(len(pelec))
    Ebatt = np.zeros(len(pelec))
    Prheos = np.zeros(len(pelec))
    VtrainBatt = np.zeros(len(pelec))
    VtrainBatt[0]=VSST 
    # EbattMAX=capacite_batterie*VSST # energie batterie = capacite * v nominale et v nominale doit etre environ égale à la tension de LAC
    Ebatt0 = EbattMAX*3/4 # comment tu as choisi ca ? énoncé ?
    Ebatt[0] = Ebatt0

    for i in range(1, len(pelec)):
        # Loi de gestion de la batterie
        if (pelec[i]<0 or V[i] == 0) and Ebatt[i-1] < EbattMAX:
            Pbatt[i] = abs(pelec[i])
            Ebatt[i] = Ebatt[i-1] + Pbatt[i]*1/3600
            if Ebatt[i] > EbattMAX:
                Ebatt[i] = EbattMAX
                Prheos[i] = Prheos[i-1] + pelec[i] + (Ebatt[i]-EbattMAX)*3600
        elif pelec[i]> seuil and Pbatt[i-1] > 0:
            Pbatt[i] = Ebatt[i-1] * 3600
            if pelec[i] > Pbatt[i]:
                pelec[i] -= Pbatt[i]
            else:
                Ebatt[i] = Ebatt[i-1] - pelec[i]*1/3600
                pelec[i] = 0
        else:
            Ebatt[i] = Ebatt[i-1]
            Pbatt[i] = Pbatt[i-1]

        PLAC[i] = pelec[i] - Pbatt[i] + Prheos[i]
        if PLAC[i] < 0:
            PLAC[i] = 0
        # print("Pbatt =", Pbatt[i], "Ebatt =", Ebatt[i], "PLAC =", PLAC[i])
        racine = VSST**2 - 4*Req[i]*(PLAC[i]/0.8)
        racine = max(racine,0)
        vtrain = (VSST + np.sqrt(racine))/2
        VtrainBatt[i] = vtrain
        # print(VtrainBatt[i])

    return VtrainBatt,Ebatt,Pbatt

Pelec=remplissage_P_elec(Pm)

def find_non_dominated_solution(objectif1, objectif2,nbre_simulations):
    solutions_non_dominees = []
    for i in range(nbre_simulations):
        is_dominated = False
        for j in range(nbre_simulations):
            if (objectif1[j] <= objectif1[i] and objectif2[j] <= objectif2[i] and 
                (objectif1[j] < objectif1[i] or objectif2[j] < objectif2[i])):
                is_dominated = True
                break
        if not is_dominated:
            solutions_non_dominees.append(i)

    return solutions_non_dominees

def monte_carlo(nbre_simulations,capacite_batterie_random,seuil_random,pelec):
    dv_max =[]
    capacite_retour = []
    vtrainbatt=np.zeros(len(pelec))
    for i in range(nbre_simulations):

        for j in range(nbre_simulations):
            vtrainbatt=np.zeros(len(pelec))
            vtrainbatt, _, _ = gestion_batterie(pelec, capacite_batterie_random[i], seuil_random[j])
            dv_max.append(VSST - np.min(vtrainbatt))
            capacite_retour.append(capacite_batterie_random[i])
            # print(dv_max[-1])
        # vtrainbatt=np.zeros(len(pelec))
    return dv_max, capacite_retour

def gestion_batterie_var(pelec,EbattMAX,seuil, v = V, plac = PLAC, vsst = VSST):
    # Initialisation
    Pbatt = np.zeros(len(pelec))
    Ebatt = np.zeros(len(pelec))
    Prheos = np.zeros(len(pelec))
    VtrainBatt = np.zeros(len(pelec))
    VtrainBatt[0]=vsst
    # EbattMAX=capacite_batterie*VSST # energie batterie = capacite * v nominale et v nominale doit etre environ égale à la tension de LAC
    Ebatt0 = EbattMAX*3/4 # comment tu as choisi ca ? énoncé ?
    Ebatt[0] = Ebatt0

    for i in range(1, len(pelec)):
        # Loi de gestion de la batterie
        if (pelec[i]<0 or v[i] == 0) and Ebatt[i-1] < EbattMAX:
            Pbatt[i] = abs(pelec[i])
            Ebatt[i] = Ebatt[i-1] + Pbatt[i]*1/3600
            if Ebatt[i] > EbattMAX:
                Ebatt[i] = EbattMAX
                Prheos[i] = Prheos[i-1] + pelec[i] + (Ebatt[i]-EbattMAX)*3600
        elif pelec[i]> seuil and Pbatt[i-1] > 0:
            Pbatt[i] = Ebatt[i-1] * 3600
            if pelec[i] > Pbatt[i]:
                pelec[i] -= Pbatt[i]
            else:
                Ebatt[i] = Ebatt[i-1] - pelec[i]*1/3600
                pelec[i] = 0
        else:
            Ebatt[i] = Ebatt[i-1]
            Pbatt[i] = Pbatt[i-1]

        plac[i] = pelec[i] - Pbatt[i] + Prheos[i]
        if plac[i] < 0:
            plac[i] = 0
        # print("Pbatt =", Pbatt[i], "Ebatt =", Ebatt[i], "PLAC =", PLAC[i])
        racine = vsst**2 - 4*Req[i]*(plac[i]/0.8)
        racine = max(racine,0)
        vtrain = (vsst + np.sqrt(racine))/2
        VtrainBatt[i] = vtrain
        # print(VtrainBatt[i])

    return VtrainBatt,Ebatt,Pbatt

# Paramètres du problème
N = 100
  # Nombre de points d'échantillonnage
E_bat_min, E_bat_max = 0, 140000  # Capacité de la batterie (en kWh)
P_seuil_min, P_seuil_max = 0, 1e6   # Puissance seuil (en MW)
V_nominal = 790  # Tension nominale du système (V)
V_min = 500  # Tension minimale admissible du train (V)

# Génération des échantillons (1000 couples capacité / puissance seuil)
np.random.seed(45)  # Pour la reproductibilité
E_bat = np.random.uniform(E_bat_min, E_bat_max, N)  # Capacités de batterie
# E_bat = np.arange(E_bat_min, E_bat_max, (E_bat_max-E_bat_min)/N)
P_seuil = np.random.uniform(P_seuil_min, P_seuil_max, N)   # Puissances seuil (en MW pour correspondre aux axes)

# Simuler les performances pour chaque configuration
P_train = np.random.uniform(0, 1e6, len(Pelec))  # Simulation de la puissance demandée (en MW)
delta_V_max = np.zeros(N)  # Chute de tension maximale pour chaque configuration
valid = np.zeros(N, dtype=bool)  # Tableau de validité des configurations


ValidDV = []
ValidCapacite = []
DV_Valid = []

for i in range(N):
    # for j in range(N):
    ValidVtrain, _, _ = gestion_batterie_var(P_train, E_bat[i], P_seuil[i])
    ValidDV = np.min(ValidVtrain)
    ValidCapacite.append(E_bat[i])
    DV_Valid.append(VSST - ValidDV)

Solutions_non_dominees=find_non_dominated_solution(E_bat , P_seuil,N)

# Affichage des solutions
plt.subplot(211)
plt.scatter(E_bat , P_seuil, color = 'skyblue')
# plt.scatter(capacite_correcte,seuil_correcte,color='red')
plt.xlabel('Capacité en énergie de la batterie (kWh)')
plt.ylabel('P seuil (MW)')
plt.title('Espace des solutions / de recherche')
plt.grid()
plt.legend()


plt.subplot(212)
plt.scatter(E_bat, DV_Valid, color = 'skyblue')
for ii in range(0, len(Solutions_non_dominees)):
    plt.scatter(ValidCapacite[Solutions_non_dominees[ii]], DV_Valid[Solutions_non_dominees[ii]], color='red', label='Solutions non dominées')
plt.xlabel('Capacité en énergie de la batterie (kWh)')
plt.ylabel('dV max (V)')
plt.title('Espace des objectifs')
plt.grid()
plt.legend()


#%%
# # Calculer la chute de tension pour chaque couple de (capacité, puissance seuil)
# for i in range(N):
#     C_batt = E_bat[i]  # Conversion de la capacité de la batterie en Wh (kWh -> Wh)

#     # Si la puissance demandée est inférieure au seuil, on n'utilise que la LAC
#     if P_train[i] <= P_seuil[i]:
#         V_train = V_nominal - (P_train[i] / C_batt)  # La chute de tension pour P_train < P_seuil
#     else:
#         # Si la puissance demandée est supérieure au seuil, on utilise la batterie et la LAC
#         V_train = V_nominal - (P_seuil[i] / C_batt)  # La chute de tension pour P_train >= P_seuil

#     # Calcul de la chute de tension maximale
#     delta_V_max[i] = V_nominal - V_train

#     # Vérifier si la configuration est valide (si la tension est suffisante)
#     if V_train >= V_min:
#         valid[i] = True

# # Extraire les performances des configurations valides
# E_bat_valid = E_bat[valid]
# delta_V_max_valid = delta_V_max[valid]
# P_seuil_valid = P_seuil[valid]

# # Identifier les solutions non dominées (Pareto) dans l'espace des objectifs
# def is_pareto_efficient(costs):
#     is_efficient = np.ones(costs.shape[0], dtype=bool)
#     for i, c in enumerate(costs):
#         if is_efficient[i]:
#             is_efficient[is_efficient] = np.any(costs[is_efficient] < c, axis=1)
#             is_efficient[i] = True
#     return is_efficient

# # Construire l'espace des objectifs
# objectives = np.vstack((E_bat_valid, delta_V_max_valid)).T
# pareto_mask_objectives = is_pareto_efficient(objectives)

# # Solutions non dominées dans l'espace des objectifs (Pareto)
# E_bat_pareto_objectives = E_bat_valid[pareto_mask_objectives]
# delta_V_max_pareto_objectives = delta_V_max_valid[pareto_mask_objectives]

# # Identifier les solutions non dominées (Pareto) dans l'espace des solutions (capacité, puissance seuil)
# solutions = np.vstack((E_bat_valid, P_seuil_valid)).T
# pareto_mask_solutions = is_pareto_efficient(solutions)

# # Solutions non dominées dans l'espace des solutions (Pareto)
# E_bat_pareto_solutions = E_bat_valid[pareto_mask_solutions]
# P_seuil_pareto_solutions = P_seuil_valid[pareto_mask_solutions]

# # Trouver la "meilleure" solution parmi les solutions Pareto dans l'espace des objectifs
# # Critère : minimiser la chute de tension et maximiser la capacité de la batterie
# best_solution_idx = np.argmin(delta_V_max_pareto_objectives)  # Choisir le meilleur selon la chute de tension

# # Meilleure solution
# best_E_bat = E_bat_pareto_objectives[best_solution_idx]
# best_delta_V_max = delta_V_max_pareto_objectives[best_solution_idx]

# # Affichage des résultats

# # 1. Affichage de l'espace de solution (capacité et puissance seuil)
# plt.figure(figsize=(15, 6))

# # Espace des solutions
# plt.subplot(1, 2, 1)
# plt.scatter(E_bat_valid, P_seuil_valid, color='blue', alpha=0.5, label="Solutions valides")
# plt.scatter(E_bat_pareto_solutions, P_seuil_pareto_solutions, color='orange', label="Solutions Pareto")
# plt.scatter(best_E_bat, P_seuil_pareto_solutions[best_solution_idx], color='green', label="Meilleure solution", s=100)
# plt.xlabel("Capacité de la batterie (kWh)")
# plt.ylabel("Puissance seuil (MW)")
# plt.title("Espace des solutions")
# plt.legend()
# plt.grid(True)

# # Espace des objectifs avec ajustement de l'axe des ordonnées pour une meilleure visualisation
# plt.subplot(1, 2, 2)
# plt.scatter(E_bat_valid, delta_V_max_valid, label="Solutions valides", alpha=0.3, color='blue')
# plt.scatter(E_bat_pareto_objectives, delta_V_max_pareto_objectives, color='orange', label="Solutions Pareto")
# plt.scatter(best_E_bat, best_delta_V_max, color='green', label="Meilleure solution", s=100)
# plt.xlabel("Capacité de la batterie (kWh)")
# plt.ylabel("Chute de tension maximale (V)")
# plt.title("Espace des objectifs")
# plt.ylim(0, 400)  # Ajustez cette valeur selon vos résultats
# plt.legend()
# plt.grid(True)

# # Affichage global
# plt.tight_layout()
# plt.show()