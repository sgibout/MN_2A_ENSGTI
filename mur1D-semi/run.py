#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mur (1D cartésien) instationnaire - Semi-implicit

@author: sgibout
"""

import numpy as np
import matplotlib.pyplot as plt

# ========================================
# Déclaration des paramètres
# ========================================

mMax = 20  # [-] nombre de pas d'espace
iMax = 1000  # [-] nombre de pas de temps

L = 0.1  # [m] épaisseur du mur
D = 360  # [s] durée de la simulation

rho = 500  # [kg/m3] masse volumique
C = 800  # [kg/m3] capacité calorifique
k = 1.0  # [W/(m.K)] conductivité thermique

TG = 30  # [°C] température imposée en x=0
TINF = 10  # [°C] température du fluide "loin" de la paroi
H = 200  # [W/(m2.K)] coefficient d'échange convectif

T0 = 20  # [°C] température initiale

# ========================================
# Construction du système matriciel
# ========================================

# Attention conernant les indices: en programmation Python/Numpy, les indices débutes à 0...
# on aura donc un décalage "de 1" sur les indices d'esapce notamment


dx = L / mMax
dt = D / (iMax - 1)  # i varie de 0 à iMax-1 et D=dt*(iMax-1)

# Ici j'introduit le Heq qui n'est autre que 1/R = 1/(1/H + 0.5*dx/k)
Heq = 1 / (1 / H + 0.5 * dx / k)

# Les constantes
alpha = k * dt / (rho * C * dx ** 2)
beta = Heq * dt / (rho * C * dx)

# Déclarations des matrices et vecteurs : forme A*T(i+1) = B*T(i)+V
A = np.zeros((mMax, mMax))
B = np.zeros((mMax, mMax))
V = np.zeros((mMax, 1))

# Construction de A - première ligne (m=1-1)
A[0, 0:2] = [1 + 1.5 * alpha, -0.5 * alpha]
# Construction de A - dernière ligne (m=mMax-1)
A[-1, -2:] = [-0.5 * alpha, 1 + 0.5 * alpha + 0.5 * beta]
# Construction de A - les autres lignes (m=2-1...mMax-1-1)
for m in range(1, mMax - 1):
    A[m, m - 1:m + 2] = [-0.5 * alpha, 1 + alpha, -0.5 * alpha]

# Construction de B (même principe)
B[0, 0:2] = [1 - 1.5 * alpha, 0.5 * alpha]
B[-1, -2:] = [0.5 * alpha, 1 - 0.5 * alpha - 0.5 * beta]
for m in range(1, mMax - 1):
    B[m, m - 1:m + 2] = [0.5 * alpha, 1 - alpha, 0.5 * alpha]

# vecteur V
V[[0, -1], 0] = [2 * alpha * TG, beta * TINF]

# ========================================
# Résolution
# ========================================

# Pour ce premier code, on va "simplement" inverser le système A*T(i+1) = B*T(i)+V
# pour calculer T(i+1) directement et on stockera TOUS les pas de temps (ce qui n'est en général pas souhaitable/possible)
# Je vais donc construire une "matrice T" qui contient les vecteurs colonnes [T(0) , T(1)], etc...

T = np.zeros((mMax, iMax))

T[:, 0] = T0  # Condition initiale

# A*T(i+1) = B*T(i)+V devient T(i+1) = inv(A)*B*T(i)+inv(A)*V = E*T(i)+F

iA = np.linalg.inv(A)
E = iA @ B
F = iA @ V

# Boucle de résolution
for i in range(0, iMax - 1):
    T[:, i + 1:i + 2] = E @ T[:, i:i + 1] + F

# ========================================
# Tracé
# ========================================

fig = plt.figure()
ax = plt.axes(projection='3d')

x = dx * (np.arange(0, mMax) - 0.5 + 1)
t = dt * np.arange(0, iMax)
X, Y = np.meshgrid(t, x)

#surf = ax.plot_surface(X, Y, T, rstride=10, cstride=20050, cmap='jet', edgecolor='black')
surf = ax.plot_surface(X, Y, T, cmap='jet', edgecolor='black')
# surf = ax.plot_surface( X,Y,T, cmap='jet', edgecolor='black')
ax.set_title('Mur 3D');
ax.set_xlabel('$t$ [s]')
ax.set_ylabel('$x$ [m]')
ax.set_zlabel('$T(x,t)$ [°C]')

fig.colorbar(surf, shrink=0.5, aspect=4)

plt.show()



# ========================================
# Petites vérifications
# ========================================

# Si on a atteint le régime permanant, alors la densité de flux doit êter constante partout

Tfinal = T[:, -1]

# A chaque interface (interne):
dens_flux = k * (Tfinal[1:-1] - Tfinal[0:-2]) / dx

R = 1 / H + L / k  # attention, pas le Heq ici
dens_flux_theo = abs(TINF - TG) / R

print(
    f"Densité de flux moyenne aux interfaces internes : {abs(np.mean(dens_flux)):7.2f} W/m2 (sigma={np.std(dens_flux):5.3f})")
print(f"Densité de flux attendue (analogie électrique)  : {dens_flux_theo:7.2f} W/m2")


