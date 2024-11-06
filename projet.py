#%% Création des données

VSST = 790 #Volts
RSST = 33 #mOhm
RHOLAC = 131e-6 # Ohm/m
RHORAIL = 18e-6 # Ohm/m

fichier = "marche_train.m"
V = 1 #km/h-1 vitesse du train à déterminer avec le fichier marche_train.m
M = 70 #tonnes masse du train
A0 = 780 #N constante forces
A1 = 6.4 #N/tonnes constante accélération
B0 = 0.0 # constante nulle ?
B1 = 0.14 #N/tonnes/(km•h-1) constante
C0 = 0.3634 #N/tonnes/(km•h-1)^2 constante inverse vitesse
C1 = 0.0

FR = (A0 + A1*M) + (B0 + B1*M)*V + (C0 + C1*M)*V**2
