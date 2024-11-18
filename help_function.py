from matplotlib import pyplot as plt
import numpy as np

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

def get_T_X(filename):
    Times = []
    X = []
    
    file = open(filename, 'r')
    readedlist = readlist(file)
    
    while readedlist:
        Time, Position = readedlist
        Times.append(Time)
        X.append(Position)
        readedlist = readlist(file)
    return np.array(Times), np.array(X)

def get_V(Times, X):
    taille = len(X)
    
    V = []
    
    V.append(0) #TODO revoir différence fini ordre supérieure
    for i in range(taille - 1):
        if i<2 or i>taille-3:
            v = (X[i+1]-X[i])/(Times[i+1]-Times[i])
        else:
            v = (X[i-2] - 8*X[i-1] + 8*X[i+1] - X[i+2])/(12*(Times[i+1]-Times[i]))
        V.append(v)

    V = np.array(V) 
    return V

def get_Acc(Times,V):
    taille = len(V)
    
    Acc = []
    
    Acc.append(0)
    for i in range(taille-1):
        if i<2 or i>taille-3:
            a = (V[i+1]-V[i])/(Times[i+1]-Times[i])
        else:
            a = (V[i-2] - 8*V[i-1] + 8*V[i+1] - V[i+2])/(12*(Times[i+1]-Times[i]))
        Acc.append(a)

    Acc = np.array(Acc)
    return Acc