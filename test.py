import numpy as np



def croisement(parent1, parent2, nombre_croisement): # le rate 0.5 signifie une chance égale , 50% des cas --> croisement réalisé
    
    
    new_individus=[] # a ajouter dans la nouvelle population apres l'appel des fonctions
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
            new.append(parent1[1])
        elif indice_capacite==3:
            new.append(parent2[1])
            
        new_individus.append(new)
      
    return new_individus


result=croisement([3,7],[2,4],10)

print(result)