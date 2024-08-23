import sys
from cProfile import label

import numpy as np
from matplotlib import pyplot as plt

from caes.air import State

# Rappel des notations
# 1 -> Extérieur / entrée compresseur
# 2 -> Sortie isentropique compressseur
# 3 -> Sortie réelle compressseur  / entrée cuve
# 4 -> Cuve / entrée turbine
# 5 -> Sortie isentropique turbine
# 6 -> Sortie réelle trubine (vers extérieur)


# Configuration

# Contient 3 colonnes (t en jour, GHI en W/m2 et Conso en W) au pas d'1 minute

#NMAX = len()
datas = np.load("data.npy")[1:24*60*7,:]

time = datas[:, 0]
GHI = datas[:, 1]
CONSO = datas[:, 2]

# Caractértistiques de la cuve
R_CUVE = 1.2 # Rayon de la cuve [m]
L_CUVE = 1 # Longeur de la partie centrale de la cuve[m]
alpha = 0.5 # coefficient d'échange global [W/m2/K]


# Conditions extérieures (supposées constantes)
Text = 20  # Température extérieure [°C]
Pext = 1  #  Pression atmosphérique [bar]

# Production solaire
S_PV = 25  # Surface de captation solaire [m2]
eta_S = 0.2  # Rendement conversion PV [-]


# Compresseur
vol_bal = 1e-2  # [m3/s]
eta_comp = 0.9 # rendement isentropique

# Turbine
eta_turb = 0.9 # rendement isentropique


# Paramètres de régulation
P_max = 80  # bar
P_min = 10
beta = 0.8 # limite pour démarrage compresseur

# ===============================
# Calculs préliminaires
# ===============================

# Volume de la cuve
V_cuve = 4*np.pi*R_CUVE**3/3 + L_CUVE*np.pi*R_CUVE**2
print(f"Vcuve = {V_cuve}m3")
# Surface de la cuve
S_cuve = 4*np.pi*R_CUVE**2 + L_CUVE*2*np.pi*R_CUVE
print(f"Scuve = {S_cuve}m3")
# ===============================
# Initialisation
# ===============================

# Calcul du dt à partir des données expérimentales
tmp = time[1:]-time[0:-1]
dt = np.mean(tmp)
if np.std(tmp) > dt/1000: # par exemple
    print("Erreur -> pas de temps non constant")
    sys.exit(-1)
dt = dt*24*60*60 # On bascule en seconde...
print(f"dt={dt} s")

# Px correspond à l'état du point x à l'instant courant [i]
P1 = State().set_T_P(Text, Pext)
P2 = State().set_T_P(Text, Pext)
P3 = State().set_T_P(Text, Pext)
P4 = State().set_T_P(Text, Pext)
P5 = State().set_T_P(Text, Pext)
P6 = State().set_T_P(Text, Pext)

# La masse d'air dans la cuve
Minit = M = V_cuve/P4.v
# L'énergie interne du système "air"
Uinit = U = M*P4.u

SOM_E = 0 # pour la vérifification du 1er principe
SOM_M = 0 # pour la vérifification du 1er principe

SOM_ACHAT = 0
SOM_VENTE = 0
SOM_SOL = 0
SOM_STOCK = 0
SOM_DESTOCK = 0
SOM_CONSO = 0


deb_comp = 0
deb_turb = 0

# On sotcke la pression et la température dans la cuve
Tcuve = np.zeros(len(time))
T6 = np.zeros(len(time))
T3 = np.zeros(len(time))

Pcuve = np.zeros(len(time))
Fachat = np.zeros(len(time))
Fvente = np.zeros(len(time))
Fcomp = np.zeros(len(time))
Fturb = np.zeros(len(time))
Fsol = np.zeros(len(time))
Fconso = np.zeros(len(time))

Dcomp = np.zeros(len(time))
Dturb = np.zeros(len(time))


for i in range(0, len(time)):

    # On mémorise les valeurs utiles
    Tcuve[i] = P4.T
    Pcuve[i] = P4.P

    # On calcule le productible du champs PV (modèle très simplifié)
    F_sol = S_PV * GHI[i] * eta_S

    # la puissance consommée
    F_conso = CONSO[i]

    # Pour pouvoir prendre les décisions, il est nécessaire de calculer les états P2/P3 et P5/P6
    # en cas de fonctionnement.
    # NB : on fait ici le choix de mettre à jour directement les états. Si la décision est prise de ne
    # pas faire fonctionner un équipement, il faudra les remettre à la bonne valeur.

    # Pour le compresseur
    tau = P4.P / P1.P  # pas besoin de passer en Pa
    T2 = (P1.T + 273.15) * np.exp(State.r * np.log(tau) / State.CP) - 273.15
    P2.set_T_P(T2, P4.P)
    # P3 -> réel
    # ON calcule h3
    h3 = P1.h + (P2.h - P1.h) / eta_comp
    P3.set_h_P(h3, P4.P)

    # Pour la turbine
    tau = P6.P / P4.P
    T5 = (P4.T + 273.15) * np.exp((State.R / State.Mmol) * np.log(tau) / State.CP) - 273.15
    P5.set_T_P(T5, P6.P)
    # P6 -> réel
    # ON calcule h6
    h6 = P4.h + eta_turb * (P5.h - P4.h)
    P6.set_h_P(h6, P6.P)

    # ON paut maintenant prendre les décisions
    # On calcule la puissance consommée par le compresseur si on le démarre
    deb_comp = vol_bal * np.exp(-0.05 * (P3.P / P1.P - 1))  # volumique
    deb_comp = deb_comp / P1.v
    F_comp = deb_comp * (P3.h - P1.h)

    F_turb = 0

    # Prise de décision
    F_vente = F_achat = 0
    F_ex = F_sol - F_conso
    if F_ex > 0:
        F_turb = 0 # on ne turbine pas car exces d'énergie
        # On a un excédent.
        if P4.P < P_max:
            # On peut comprimer
            if F_ex > beta * F_comp:
                # On démarre le compresseur
                # On garde F_comp et deb_comp
                F_achat = max(0,F_comp-F_ex)
                F_vente = max(0,F_ex-F_comp)
            else:
                # On ne démarre pas le compresseur et on vend tout
                F_comp = 0
                deb_comp = 0
                F_vente = F_ex
        else:
            # On ne démarre pas le compresseur et on vend tout
            deb_comp = 0
            F_comp = 0
            F_vente = F_ex
    else:
        # On a un déficit... Détente éventuelle
        deb_comp = 0
        F_comp = 0
        if P4.P > P_min:
            # ON peut détendre
            F_turb = np.abs(F_ex)
        else:
            # On ne peut pas lancer la turbine et on doit acheter
            F_turb = 0
            F_achat = np.abs(F_ex)

    # Calcul de la puissance et du débit massique turbine
    if P4.h>P6.h:
        deb_turb = F_turb / (P4.h - P6.h)
    else:
        deb_turb = 0

    # On connait mantenant le mode de fonctionnement
    # notamment les puissances et les débits entrée/sortie compresseur et turbine

    #print(f"{time[i]}\n  P1={P1}\n  P2={P2}\n  P3={P3}\n  P4={P4}\n  P5={P5}\n  P6={P6}")

    # On calcule les états autres que P4
    # P1 ne change pas (sauf si on introduit des conditions météos variables

    # P2 et P3 sont égaux à P4 si le compresseur est à l'arret. Sinon on suit l'algo du poly
    if deb_comp==0:
        # C'est ici qu'on <<corrige>> P2/P3 si le compresseur ne fonctionne pas
        P2.copy(P4)
        P3.copy(P4)

    # P5 et P6 se calcule avec la détente...
    if deb_turb == 0:
        # Ici le choix est plus difficile. Quel l'état à la sortie de la turbine lorsqu'elle ne
        # fonctionne pas ? On peut apr exemple considérer que la canalisation se remplie d'air
        # extérieur. Ce choix n'a cependant aucune influence sur le déroulement de la simualtion
        # (regarder vos équations)
        P5.copy(P1)
        P6.copy(P1)

    # On fait évoluer l'état de l'air dans la cuve -> bilans
    F_ech = alpha*S_cuve*(Text-P4.T)
    # Bilan masse
    delta_mass = dt*(deb_comp - deb_turb)
    # Bilan énergie
    delta_energy = dt * (F_ech + deb_comp * P3.h - deb_turb * P4.h)

    SOM_E += delta_energy
    SOM_M += delta_mass

    SOM_ACHAT += dt*F_achat
    SOM_VENTE += dt*F_vente
    SOM_STOCK += dt*F_comp
    SOM_DESTOCK += dt*F_turb
    SOM_SOL += dt * F_sol
    SOM_CONSO += dt * F_conso

    # On calcule le nouvel état (à la fin du pas de temps actuel)
    M = M + delta_mass
    U = U + delta_energy

    # P4 se calcule à partir de M et U
    u = U/M # Energie interne massique -> nous donne
    v = V_cuve/M
    P4.set_u_v(u,v)

    T6[i] = P6.T
    T3[i] = P3.T
    Fachat[i] = F_achat
    Fvente[i] = F_vente
    Fcomp[i] = F_comp
    Fturb[i] = F_turb
    Fsol[i] = F_sol
    Fconso[i] = F_conso
    Dcomp[i] = deb_comp
    Dturb[i] = deb_turb


Ufinal = M*P4.u
Mfinal = M

print(f"Energie : Ufinal-Uinit={Ufinal-Uinit} ;  SOM_E={SOM_E} ; Delta = {Ufinal-Uinit-SOM_E}")
print(f"Masse : Mfinal-Minit={Mfinal-Minit} ;  SOM_M={SOM_M} ; Delta = {Mfinal-Minit-SOM_M}")

# Affichage Bilan
print("======= BILAN ==================================")
print(f"CONSO   = {SOM_CONSO/3600/1000:.0f} kWh")
print(f"ACHAT   = {SOM_ACHAT/3600/1000:.0f} kWh")
print(f"VENTE   = {SOM_VENTE/3600/1000:.0f} kWh")
print(f"STOCK   = {SOM_STOCK/3600/1000:.0f} kWh")
print(f"DESTOCK = {SOM_DESTOCK/3600/1000:.0f} kWh")
print(f"SOLAIRE = {SOM_SOL/3600/1000:.0f} kWh")
print("------------------------------------------------")
print(f"Rendement stock    = {100*SOM_DESTOCK/SOM_STOCK:.1f}%")
print(f"Couverture solaire = {100*SOM_SOL/SOM_CONSO:.1f}%")
print(f"Taux vente/Solaire = {100*SOM_VENTE/SOM_SOL:.1f}%")
print(f"Taux achat/conso   = {100*SOM_ACHAT/SOM_CONSO:.1f}%")

print("================================================")


fig, (axis_T, axis_P, axis_conso, axis_meca, axis_flux_ext) = plt.subplots(5,1)

axis_P.plot(time,Pcuve, label="$P_4$")
axis_P.set_xticklabels([])
axis_P.set_ylabel('(Bar)')
axis_P.legend(loc='center left', bbox_to_anchor=(1, 0.5))

axis_T.plot(time,Tcuve, label='$T_4$')
axis_T.plot(time,T6, ":", label='$T_6$')
axis_T.plot(time,T3,":", label='$T_3$')
axis_T.set_xticklabels([])
axis_T.set_ylabel('(°C)')
axis_T.legend(loc='center left', bbox_to_anchor=(1, 0.5))

axis_conso.plot(time,-Fconso/1000, label="$\mathcal{P}_{conso}$")
axis_conso.plot(time,Fsol/1000, label = "$\mathcal{P}_{sol}$")
axis_conso.legend(loc='center left', bbox_to_anchor=(1, 0.5))
axis_conso.set_xticklabels([])
axis_conso.set_ylabel('(kW)')

axis_meca.plot(time,-Fcomp/1000, label="$\mathcal{P}_{comp}$")
axis_meca.plot(time,Fturb/1000, label="$\mathcal{P}_{turb}$")
axis_meca.legend(loc='center left', bbox_to_anchor=(1, 0.5))
axis_meca.set_xticklabels([])
axis_meca.set_ylabel('(kW)')


axis_flux_ext.plot(time,Fvente/1000, label="$\mathcal{P}_{vente}$")
axis_flux_ext.plot(time,-Fachat/1000, label="$\mathcal{P}_{achat}$")
axis_flux_ext.legend(loc='center left', bbox_to_anchor=(1, 0.5))
axis_flux_ext.set_xlabel('Temps (Jour)')
axis_flux_ext.set_ylabel('(kW)')


#axis_flux_ext.plot(time,Fvente, label("Vente"))
#axis_flux_ext.set_xlabel('Time (Jour)')
#axis_flux_ext.set_ylabel('(W)')
#axis_flux_ext.legend()

#plt.subplots_adjust(hspace=0.4)
plt.tight_layout(rect=[0, 0, 0.95, 1])
plt.show()









