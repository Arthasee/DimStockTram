"""Fichier pour le projet d'Optimisation
"""
#%% Imports

import random
import copy
import numpy as np
import matplotlib.pyplot as plt
from help_function import *

#%% Création et conversion des données
NB_SIMU= 1000
NB_SIMU= 1000
VLAC = 750 #Volts
VSST = 790 #Volts
RSST = 33*1e-3 #mOhm
RHOLAC = 95e-6#131e-6 # Ohm/m
RHORAIL = 10e-6#18e-6 # Ohm/m

Times = []
X = []
V = [] #km/h-1 vitesse du train à déterminer avec le fichier marche_train.txt
Acc = [] #km/h-2 accélération du train à déterminer avec le fichier marche_train.txt
RLAC1 = [] # Résistance de la LAC entre la sous-station 1 et le train (valeurs dépendante de x)
RLAC2 = [] # Résistance de la LAC entre la sous-station 2 et le train (valeurs dépendante de x)
Rrail1 = [] # Résistance du rail entre la sous-station 1 et le train (valeurs dépendante de x)
Rrail2 = [] # Résistance du rail entre la sous-station 2 et le train (valeurs dépendante de x)
# Résistance équivalente pour la partie supérieure du schéma de Thévenin,
# traversée par le courant I1 (dépend de x)
R1 = []
# Résistance équivalente pour la partie inférieure du schéma de Thévenin,
# traversée par le courant I2(dépend de x)
R2 = []
Req = [] # Résistance équivalente totale du schéma de Thévenin (dépend de x)
PLAC = [] #Puissance de la LAC (dépend de x)
Vtrain = [] #Tension du train à tout instant (dépend de x)
Itrain = [] #Intensité aux bornes du train à tout moment (dépend de x)
I1 = [] #Intensité de la partie supérieure du schéma de Thévenin
I2 = [] #Intensité de la partie inférieure du schéma de Thévenin

ALPHA = 0 # angle de la pente du chemin
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

Fm = M*Acc + M*9.81*np.sin(ALPHA) + FR # Force mécanique - ici alpha = 0

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
# On sait d'après la loi des mailles que V1 - V2 = 0, donc V1 = V2,
# donc R1 * I1 = R2 * I2, donc I2 = (R1 * I1)/R2

# D'après la loi des noeuds, I1 + I2 = Itrain, donc en remplaçant I2
# par son expression en fonction de I1 on obtient :

# I1 + (R1 * I1)/R2 = Itrain, donc I1(R2 + R1)/R2 = Itrain,
# donc I1 = (R2 * Itrain)/(R1 + R2)
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

def remplissage_p_elec(pm):
    """Fonction permettant de remplir la puissance électronique nécessaire

    Args:
        pm (Array): Puissance mécanique

    Returns:
        Array: Puissance électronique
    """

    # Initialisation
    pelec = np.zeros(len(pm))

    # Remplissage de Pelec, rendement de perte 0.8
    for i in range(1, len(pm)):
        if pm[i]<=0:
            pelec[i] = pm[i]*0.8
        else:
            pelec[i] = pm[i]/0.8
    return pelec

def gestion_batterie(pelec,ebatt_max,seuil, req = Req, vsst = VSST, plac = PLAC):
    """Fonction calculant la gestion de la batterie pour chaque pas de temps

    Args:
        pelec (Array): Puissance électronique nécessaire
        EbattMAX (int): Capacité de la Batterie
        seuil (int): Seuil où la batterie doit intervenir
        v (Array, optional): Vitesse. Defaults to V.
        vsst (int, optional): Tension nominale. Defaults to VSST.
        plac (Array, optional): Puissance de la LAC. Defaults to PLAC.

    Returns:
        Tuple: Liste de la Tension du Train, de l'Energie de la batterie et de sa puissance
    """
    # Initialisation
    pbatt = np.zeros(len(pelec))
    ebatt = np.zeros(len(pelec))
    p_rheos = np.zeros(len(pelec))
    v_train_batt = []
    v_train_batt.append(vsst)
    ebatt0 = ebatt_max*4/4 # comment tu as choisi ca ? énoncé ?
    ebatt[0] = ebatt0
    pbatt[0] = ebatt[0]*3600
    compt_if = 0
    compt_elif = 0
    compt_else = 0

    for i in range(1, len(pelec)):
        # Loi de gestion de la batterie
        if (pelec[i]<0 ) and ebatt[i-1] <= ebatt_max: #or v[i] == 0
            compt_if += 1
            pbatt[i] = abs(pelec[i])
            ebatt[i] = ebatt[i-1] + pbatt[i]*1/3600
            if ebatt[i] > ebatt_max:
                p_rheos[i] = p_rheos[i-1] + pelec[i] + (ebatt[i]-ebatt_max)*3600
                ebatt[i] = ebatt_max
        elif pelec[i]> seuil and pbatt[i-1] > 0:
            compt_elif +=1
            pbatt[i] = ebatt[i-1] * 3600
            if pelec[i] > pbatt[i]:
                pelec[i] -= pbatt[i]
            else:
                ebatt[i] = ebatt[i-1] - pelec[i]*1/3600
                pelec[i] = 0
        else:
            compt_else += 1
            ebatt[i] = ebatt[i-1]
            pbatt[i] = pbatt[i-1]

        plac[i] = pelec[i] - pbatt[i] + p_rheos[i]
        if plac[i] < 0:
            plac[i] = 0
        # print("Pbatt =", Pbatt[i], "Ebatt =", Ebatt[i], "PLAC =", PLAC[i])
        racine = vsst**2 - 4*req[i]*(plac[i]/0.8)
        # print(plac[i])
        # print(req[i])
        racine = max(racine,0)
        vtrain = (vsst + np.sqrt(racine))/2
        v_train_batt.append(vtrain)
    # print(v_train_batt)
        # print(VtrainBatt[i])
    # print(f'if {compt_if}, elif {compt_elif}, else {compt_else}, {np.min(v_train_batt)}')
    return v_train_batt,ebatt,pbatt

Pelec=remplissage_p_elec(Pm)
EbattMax = 18*1e3

Seuil = 0.5 * np.max(Pelec) # choisi par Baptiste, quel est l'unité de Pelec à discuter
VTrainBatt,EBatt,PBatt=gestion_batterie(Pelec,EbattMax,Seuil)


# Affichage des solutions
# trace(Times, EBatt, "Temps[s]", "Energie de la batterie",
#        "Energie de la batterie en fonction du temps")
# trace(Times, PLAC, "Temps[s]", "PLAC", "PLAC avec batterie en fonction du temps")
# trace(Times, PBatt, "Temps[s]", "puissance batterie", "puissance batterie en fonction du temps")
# trace(Times, VTrainBatt, "Temps[s]", "Vtrain",
#         "Vtrain avec batterie en fonction du temps") #, [0, 140]


#%% Dimmensionnement du système de stockage
# construction de l’ensemble des solutions non dominées pour les critères
# « Capacité en énergie de la batterie » et « Chute de tension maximale »
# (qui est la différence entre la tension nominale (750V) et
# la tension réelle mesurée aux bornes du train.)

#%% Méthode de Monte-Carlo

# Construire les solutions non dominées

def find_non_dominated_solution(objectif1, objectif2,nbre_simulations):
    """Fonction trouvant les solutions non-dominées

    Args:
        objectif1 (Array): Premier objectif à analyser
        objectif2 (Array): Deuxième objectif à analyser
        nbre_simulations (int): Nombre de simulation

    Returns:
        Array: Ensemble des solutions non-dominées
    """
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

# # Stockage des résultats

# # Capacité de la batterie (en kWh) objectif1
# capacite_batterie = np.random.uniform(50, 200, NB_SIMU)
# chute_tension = np.random.uniform(10, 250, NB_SIMU)  # Chute de tension maximale (en V) objectif2
# # Appel de la fonction
# solutions_non_dominees=find_non_dominated_solution(capacite_batterie ,chute_tension,NB_SIMU)

# # Affichage des solutions
# plt.scatter(capacite_batterie, chute_tension, color = 'skyblue',
#               label='Ensemble des solutions par la méthode de Monté - Carlo')
# plt.scatter(capacite_batterie[solutions_non_dominees],
#               chute_tension[solutions_non_dominees], color='red', label='Solutions non dominées')

# # comment ca marche? plot une list de list
# plt.xlabel('Capacité en énergie de la batterie (kWh)')
# plt.ylabel('Chute de tension maximale (V)')
# plt.legend()
# plt.show()


#%% Méthode de Monte-Carlo 2

def monte_carlo(nbre_simulations,capacite_batterie_random,seuil_random, pelec, vsst = VSST, placc = PLAC, vlac = VLAC, rlac1 = RLAC1, rlac2 = RLAC2):
    """Fonction de Monte-Carlo

    Args:
        nbre_simulations (int): Le nombre de simulation
        capacite_batterie_random (Array): Liste des différentes tailles de Capacité de la batterie
        seuil_random (Array): Liste des différents seuils de quand la batterie s'active
        pelec (Araay): Puissance électronique demandé en chaque point

    Returns:
        Tuple: La liste de la chute de tension max pour chaque essai et la liste des capacités
    """
    dv_max =[]
    for i in range(nbre_simulations):
        pelecc = remplissage_p_elec(Pm)
        placc = vlac**2/(rlac1+rlac2)
        vtrainbatt, _, _ = gestion_batterie(pelecc, capacite_batterie_random[i], seuil_random[i], plac = placc)
        dv_max.append(vsst - min(vtrainbatt))
    return dv_max

# Paramètres à optimiser sont le cout et la chute de tension dv max -->
# le cout est proportionel à la capacité,
# plus la Pseuil est petit plus la batterie rentre en compte lorsquil ne faut pas
# et plus Pseuil est grand plus il y a une chute de tension -->
# parametre a optimiser capacité et chute de tension

# Capacité de la batterie (en kWh) objectif1
Capacite_batterie_random=  np.random.uniform(0, 200000, NB_SIMU)
# Capacite_batterie_random=  np.random.uniform(0, 200000, NB_SIMU)

# Chute de tension maximale (en MW) objectif2
Seuil_random = np.random.uniform(0, 1e6, NB_SIMU)


# dV_max =[]
# vtrainbatt=np.zeros(len(Pelec))
# for i in range(NB_SIMU):
#     pelecc = remplissage_p_elec(Pm)
#     placc = VLAC**2/(RLAC1+RLAC2)
#     vtrainbatt=np.zeros(len(Pelec))
#     vtrainbatt, ebattt, pbattt = gestion_batterie(pelecc, Capacite_batterie_random[i], Seuil_random[i], plac = placc)
#     print(vtrainbatt)
#     dV_max.append(VSST - np.min(vtrainbatt))

# dV_max=monte_carlo(NB_SIMU,Capacite_batterie_random,Seuil_random,Pelec)
# Solutions_non_dominees=find_non_dominated_solution(Capacite_batterie_random ,dV_max,NB_SIMU)


# Affichage des solutions
# plt.subplot(211)
# plt.scatter(Capacite_batterie_random, Seuil_random, color = 'skyblue')
# # plt.scatter(capacite_correcte,seuil_correcte,color='red')
# for ii, sol in enumerate(Solutions_non_dominees):
#     plt.scatter(Capacite_batterie_random[sol],
#                 Seuil_random[sol], color='red')
# plt.xlabel('Capacité en énergie de la batterie (kWh)')
# plt.ylabel('P seuil (MW)')
# plt.title('Espace des solutions / de recherche')
# plt.grid()
# plt.legend()


# plt.subplot(212)
# plt.scatter(Capacite_batterie_random, dV_max, color = 'skyblue')
# for ii, sol in enumerate(Solutions_non_dominees):
#     plt.scatter(Capacite_batterie_random[sol],
#                 dV_max[sol], color='red')
# plt.xlabel('Capacité en énergie de la batterie (kWh)')
# plt.ylabel('dV max (V)')
# plt.title('Espace des objectifs')
# plt.grid()
# plt.legend()


# plt.tight_layout()
# plt.show()

joulou = np.random.randint(0, len(Capacite_batterie_random), 1)

EbattMax = Capacite_batterie_random[joulou[0]]
Seuil = Seuil_random[joulou] #np.random.randint(0, len(Seuil_random), 1)

pelecc = remplissage_p_elec(Pm)
placc = VLAC**2/(RLAC1+RLAC2)
VTrainBatt,EBatt,PBatt=gestion_batterie(pelecc,EbattMax,Seuil, plac = placc)

# trace(Times, EBatt, "Temps[s]", "Energie de la batterie", "Energie de la batterie en fonction du temps")
# trace(Times, PLAC, "Temps[s]", "PLAC", "PLAC avec batterie en fonction du temps")
# trace(Times, PBatt, "Temps[s]", "puissance batterie", "puissance batterie en fonction du temps")
# trace(Times, VTrainBatt, "Temps[s]", "Vtrain", "Vtrain avec batterie en fonction du temps") #, [0, 140]


#%% Méthode NGSA2 (non-sorted algorithm system)

""" Notes: Objectif: trouver un enssemble de solutions non-dominées appelées Paréto optimales. Identification de front de paréto 
-> ensemble de solutions où aucune solution ne peut être améliorée dans un objectif sans détériorer un autre objectif.
1) Population N choisis aléatoirement
2) Individu --> gène , enssemble des gène de l'inidividu --> chromosome
3) À chaque génération (nombre de génération choisi préalablement) : 
- Tri des individus, premier front contient les solutions non dominées( Paréto optimales), le second front contient les solutions dominées uniquement par celles du premier front, et ainsi de suite.
50 pourcent des meilleurs individus seront choisis (dans chaque front de paréto j'imagine)
- Mesure de diversité, distance de regroupement
- Opération génétique classique vu en IA TP4: Mutations et Croisement pour éxplorer d'autres espaces de recherches
4) Arret à la fin du cycle de génération


3 caractéristiques suivantes: 
-Principe de l'éllitsime
-Favorise les solutions non dominées
- Donc, dévloppe une grande variété de solution


Ref:
https://www.mechanics-industry.org/articles/meca/pdf/2010/03/mi100066.pdf
https://moodle-sciences-24.sorbonne-universite.fr/pluginfile.php/227541/mod_resource/content/1/142199_BOUKIR%20_2023_archivage.pdf THESE
https://moodle-sciences-24.sorbonne-universite.fr/pluginfile.php/225975/mod_resource/content/1/OPTIM_2425_Presentation_projet.pdf SLIDE
"""

""" Notes de la thèse: 
Explorer les zones qui paraissent prometteuses sans être bloquées
par un optimum local."""


def NGSA2(nb_generation,pop_size):
    # Initialisation
  

    mutation_rate=0.5
    
    list_capacite_batterie=  np.random.uniform(0, 200000, pop_size)
    list_chute_tension=  np.random.uniform(0, 1e6, pop_size)
    
    liste_generation = []
    population=[]
    for i in range(0,pop_size):
        population.append([list_capacite_batterie[i],list_chute_tension[i]])
    
    liste_generation.append(population)
    
    for i in range(nb_generation):
        capacite_bat=[]
        chute_ten=[]
        # print(population[0])
        # séparer la capcite et la chute de tension dans deux list distintcs pour le monte carlo et rang
        for individu in population:
            capacite_bat.append(individu[0])
            chute_ten.append(individu[1])
        
        dv = monte_carlo(len(capacite_bat), capacite_bat, chute_ten, Pelec) 
        rang_li, cap_li, seuil_li, dv_li = rang(capacite_bat, chute_ten, dv)
        
        best_list= selection(rang_li,cap_li,seuil_li,pop_size)
        nb=len(best_list[0])
        nombre_enfant_souhaite=int(nb/2) # donc moitie d'enfant creer par croisement, # possible de modifier
        new_enfant_croise=[]
        for i in range(0, nombre_enfant_souhaite):
            i_parent1=int(np.random.randint(0,len(best_list)))
            i_parent2=int(np.random.randint(0,len(best_list)))
            
            new_enfant_croise= croisement(best_list[i_parent1],best_list[i_parent2],nombre_enfant_souhaite)
        nombre_mutation_souhaite=int(nb/2) # possible de modifier
        list_a_muter = choix(best_list, nombre_mutation_souhaite)
        enfants_mute=mutation(list_a_muter,mutation_rate)
        nouvelle_gen=best_list+new_enfant_croise+enfants_mute
        
        while [] in nouvelle_gen:
            nouvelle_gen.remove([])
        population=copy.deepcopy(nouvelle_gen)
        liste_generation.append(population)
     
    return population, liste_generation



# Fonctions qu'on a besoin pour réaliser l'algorithme génétique

def choix(liste_a_choix, nombre_mutation):
    id = np.random.uniform(0, len(liste_a_choix), nombre_mutation)
    liste_a_muter = []
    for i, choixe in enumerate(liste_a_choix):
        if i in id:
            liste_a_muter.append([choixe[0],choixe[1]])
            # liste_a_muter[1].append(choixe[1])
    return liste_a_muter

def rang(capacite_batterie, chute_tension, dv_max):
    Cap_batt = []
    Cap_batt.append(capacite_batterie.copy())
    Chu_tension = []
    Chu_tension.append(chute_tension.copy())
    dv = []
    dv.append(dv_max.copy())
    rang_list = []

    while len(Cap_batt[-1]) > 0:
        test = find_non_dominated_solution(Cap_batt[-1], dv[-1],len(Cap_batt[-1]))
        rang_list.append(test)
        Cap_batt.append(np.array(np.delete(Cap_batt[-1], rang_list[-1])))
        Chu_tension.append(np.array(np.delete(Chu_tension[-1], rang_list[-1])))
        dv.append(np.array(np.delete(dv[-1], rang_list[-1])))

    return rang_list, Cap_batt, Chu_tension, dv

# rang_test, Test_Capacite, _, _ =rang(Capacite_batterie_random, Seuil_random, dV_max)

def mutation(population, mutation_rate):
    mu_rate1=mutation_rate
    mu_rate2=mutation_rate
    population_mutee=[[]]
    # pop_ret
    for pop in population:
        if np.random.random() < mu_rate1+1: # pour le sueil
            pop[0] = np.random.uniform(0, 1e6) # mutations aleatoire
            
        if np.random.random()<mu_rate2: # pour la capacite
            pop[1] = np.random.uniform(0, 200000) # mutations aleatoire
            
        population_mutee.append(pop)
    return population_mutee
  

def croisement(parent1, parent2, nombre_croisement): # le rate 0.5 signifie une chance égale , 50% des cas --> croisement réalisé
    
    
    new_individus=[[]] # a ajouter dans la nouvelle population apres l'appel des fonctions
    for i in range(0,nombre_croisement):
        # les parents 1 et 2 sont choisi aléatoirement dans la fonction principale
        new=[]
        indice_seuil= np.random.randint(0, 2) # 0 ou 1  parent 1 ou parent2
        indice_capacite= np.random.randint(2, 4) # 0 ou 1
        
        if indice_seuil==0:
            new.append(parent1[0])
        elif indice_seuil==1:
            new.append(parent2[0])
        if indice_capacite==2:
            new.append(parent1[0])
        if indice_capacite==3:
            new.append(parent2[0])

        new_individus.append(new)
      
    return new_individus




def selection(rang_l, cap_l, seuil_l, pop_size, distances = 1):
    # il faut selectionner 50% des meilleurs d'après le slide du projet
    selected_cap=[]
    selected_seuil = []
    for i, cp in enumerate(cap_l):
        cp = np.array(cp)
        seuil_l[i] = np.array(seuil_l[i])
    N=int(pop_size*0.5)  # nombre_a_selectionne
    for i, rng in enumerate(rang_l):
        if len(selected_cap) + len(rng) < N:
            for r in rng:
                selected_cap = selected_cap + [cap_l[i][r]]
                selected_seuil = selected_seuil + [seuil_l[i][r]]
        else:
            reste= N-len(selected_cap)
            j = 0
            while reste>0:
                selected_cap = selected_cap + cap_l[i][rng[j]]
                selected_seuil = selected_seuil + seuil_l[i][rng[j]]
                reste = reste - 1
            break
    selected = []
    for i, sel in enumerate(selected_cap):
        selected.append([sel,selected_seuil[i]])
         
    return selected


#test:  état selection FONCTIONNE
# selectionne_test=selection([[1,2,3,4],[5,4,9,10]],0,10)
# print(selectionne_test)




poptest_final, generation=NGSA2(10,100)
cap_ngsa2 = []
seuil_ngsa2 = []
dv_ngsa2 = []

for gen in generation:
    cap = []
    seuil = []
    dv = []
    for i in range(len(gen)):
        cap.append(gen[i][0])
        seuil.append(gen[i][1])
    cap_ngsa2.append(cap)
    seuil_ngsa2.append(seuil)
    dv_ngsa2.append(monte_carlo(len(cap),cap,seuil, Pelec))

# print(len(cap_ngsa2[0]))

plt.figure()
for i in range(len(cap_ngsa2)):
    plt.scatter(cap_ngsa2[i], seuil_ngsa2[i], label = f'Génération {i}')
plt.xlabel("Capacité")
plt.ylabel("Seuil")
plt.legend()
plt.show()

plt.figure()
for i in range(len(cap_ngsa2)):
    plt.scatter(cap_ngsa2[i], dv_ngsa2[i], label = f'Génération {i}')
plt.xlabel("Capacité")
plt.ylabel("Chute de Tension")
plt.legend()
plt.show()








# def distance_encombrement(population):
    
#     return 1


population_test = NGSA2(7, 100)
print(len(population_test[-1]))


# def voir_convergence():
    
#     return 1



""" Parametre à changer pour obtenir une meilleur convergence des front de paréto:
rate mutation, nb de mutation sur un genome, rate croisement, pop_size,  a rajoutez si vous pensez à quelque chose"""
