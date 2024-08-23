import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------
# Paramétrage (utilisateur)
# --------------------------------------------

# Géométrie de la cuve cylindrique
H = 1.5         # hauteur de la cuve [m]
R = 0.3         # rayon de la cuve [m]
r = 0.01        # rayon de l'orifice [m]

# Propriété du fluide
rho = 1000      # masse volumique du fluide [kg/m3]

# Conditions aux limites
debit_in = 1    # débit massique d'alimentation [kg/s]

# Condition initiale
h0 = 0.5*H      # hauteur initiale [m]

# Durée de la simulation
D = 2000         # durée de la simulation [s]

# Paramètres de discrétisation
iMax = 1000     # nombre de pas de temps [-]

# Paramètres additionels
g = 10          # acceleration de la pesanteur [m/s2]

# --------------------------------------------
# Déclaration des variables utiles
# --------------------------------------------

# Pas de temp [s]
dt = D/(iMax-1)

# hauteur dans la cuve en fonction du temps
h = np.zeros(iMax)

# --------------------------------------------
# Condition initiale
# --------------------------------------------
h[0] = h0
var2 = 0
# --------------------------------------------
# Résolution
# --------------------------------------------
for i in range(0,iMax-1):

    # Calcul du débit sortant
    debit_out = - rho*np.sqrt(2*g*h[i])*np.pi*r**2

    # Bilan (et passage à l'itération suivante)
    h[i+1] = max(0,h[i] + dt*(debit_in + debit_out)/(rho*np.pi*R**2))
    var2 += dt*(debit_in + debit_out)

# --------------------------------------------
# Vérification conservatif
# --------------------------------------------
var1 = rho*np.pi*R**2*((h[-1]-h[0]))
print(f'var1 = {var1}')
print(f'var2 = {var2}')


# --------------------------------------------
# Affichage
# --------------------------------------------
hmax = debit_in**2/(2*rho**2*g*np.pi**2*r**4)
print(f"hmax = {hmax}")

debit_max = rho*np.sqrt(2*g*H)*np.pi*r**2
print(f"debit_max = {debit_max}")

print(100*(var1-var2)/var2)
t = np.linspace(0,D,iMax)

plt.plot(t,h)
plt.xlabel('t [s]')
plt.ylabel('h [m]')
plt.show()

