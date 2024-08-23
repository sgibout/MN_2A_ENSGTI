import numpy as np
import matplotlib.pyplot as plt

from ecslib import charge

# --------------------------------------------
# Paramétrage (utilisateur)
# --------------------------------------------

# Géométrie de la cuve cylindrique
H = 1.5         # hauteur de la cuve [m]
R = 0.3         # rayon de la cuve [m]

# Propriété du fluide
rho = 1000      # masse volumique du fluide [kg/m3]
C = 1000        # [J/kg/K]

# Propriétés d'échange
K = 0.1         # coefficient d'échange global [W/m2/K]
Tinf = 20       # Température à l'extérieur [°C]

# Condition initiale
T0 = 20         # [°C]

Tf = 15         # Température eau froide [°C]

# Résistance électrique
TC = 60         # [°C] Consigne
DTC = 5         # [°C]
PMAX = 2000     # [W] puissance maximale (on)

# On va caler le pas et la durée de simulation sur le profil de puisage
dt, iMax , profil = charge("profil.csv")
Tp  = 45        # Température de puisage [°C]

# --------------------------------------------
# Déclaration des variables utiles
# --------------------------------------------

# Etat de fonctionnement de la résitance
on = False # à l'arret au démarrage

# Volume
V = np.pi*R**2*H
S = 2*np.pi*R*H + 2*(np.pi*R**2)

# Température
T = np.zeros(iMax)


# On mémorise pour pouvoir analyser
Pa = np.zeros(iMax) # Puissance d'appoint
Pe = np.zeros(iMax) # Puissance élecrique
dm = np.zeros(iMax) # débit extrait

# --------------------------------------------
# Condition initiale
# --------------------------------------------
T[0] = T0

v2 = 0
# --------------------------------------------
# Résolution
# --------------------------------------------
for i in range(0,iMax-1):

    # Puissance électrique
    # --------------------
    if on and T[i]>TC+DTC:
        on = False
    if not on and T[i]<TC-DTC:
        on = True

    Pe[i] = PMAX if on else 0

    # Puisance "pedue"
    # ----------------
    Pp = K*S*(Tinf-T[i])

    # Advection
    # ---------

    # on récupère le débit puisé
    dmp = profil[i]
    if T[i]<Tp:
        # on doit chauffer
        Pa[i] = dmp*C*(Tp-T[i])
        dm[i] = dmp
    else:
        # On doit refroidir
        dm[i] = dmp*(Tp-Tf)/(T[i]-Tf)

    Padv = dm[i]*C*(Tf-T[i])

    # Bilan (et passage à l'itération suivante)
    T[i+1] = T[i] + dt*(Pe[i]+Pp+Padv)/(rho*C*V)

    v2 += dt*(Pe[i]+Pp+Padv)

# --------------------------------------------
# Vérification des bilans
# --------------------------------------------
v1 = rho*C*V*(T[-1]-T[0])
print(f"Vérification du bilan : dE/E = {100*(v1-v2)/v2:.3e} % ?")


# --------------------------------------------
# Affichage
# --------------------------------------------

t = dt*np.linspace(0,iMax-1,iMax)


fig, axes = plt.subplots(3, 1)
axes[0].plot(t,T)
axes[0].plot([0, t[-1]],[Tp, Tp],'r:')
axes[0].set_ylabel('T [°C]')
axes[0].set_xticklabels([])

axes[1].plot(t,Pa, label="$P_a$")
axes[1].plot(t,Pe, label="$P_e$")
axes[1].set_ylabel('$P$ [W]')
axes[1].legend(loc='upper left')
axes[1].set_xticklabels([])

axes[2].plot(t,dm,label=r"$\dot{m}$")
axes[2].plot(t,profil,label=r"$\dot{m}_p$")
axes[2].plot(t,profil-dm,label=r"$\dot{m}_f$")
axes[2].legend(loc='upper left')
axes[2].set_ylabel(r'$\dot{m}$ [kg/s]')

plt.show()

print (3600*PMAX/rho/C/V)