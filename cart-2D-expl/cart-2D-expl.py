#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: sgibout
"""
import os

import numpy as np
import matplotlib.pyplot as plt

# ========================================
# Déclaration des paramètres
# ========================================

mMax = 30       # [-] nombre de pas d'espace selon X
nMax = 30       # [-] nombre de pas d'espace selon Y

LX = 0.1         # [m] longueur mur selon X
LY = 0.1         # [m] longueur mur selon Y
D = 3600        # [s] durée de la simulation

dt = 1       # [s] pas de temps de simulation
dtLog = 60       # [s] pas de temps de sauvegardeè2

rho = 1000       # [kg/m3] masse volumique

K = lambda T : 1 + 0.1 * (T - 20)
C = lambda T : 2000 + 20 * (T - 20)

TINF = lambda t : 20 +  10*np.sin(2*np.pi*t/1800) # [°C] température externe

H1 = 100         # [W/m2/K] Coefficient d'échange global (haut)
H2 = 100         # [W/m2/K] Coefficient d'échange global (droite)

PHI = 1000     # [W/m2] densité de flux imposée

# ========================================
# préparation des géométries
# ========================================

dx = LX/mMax
dy = LY/nMax
Sx = 1*dy
Sy = 1*dx
V = 1*dx*dy

# ========================================
# Routine d'export
# ========================================

# Création du répertoire si nécessaire
if not os.path.exists("res"):
    os.makedirs("res")

# Pour gagner un peu de temps...
X = [ (m+0.5)*dx for m in range(0,mMax) ]
Y = [ (n+0.5)*dy for n in range(0,nMax) ]
XX, YY = np.meshgrid(X,Y)

def export():
    levels = np.linspace(10, 35, num=20)
    c1=plt.contourf(XX,YY,T,     cmap ="jet", levels=levels)
    #c1=plt.pcolormesh(T,     cmap ="jet", vmin = 10, vmax = 35)
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title(f"T(x,y,t) @ t={time:.0f}s")
    plt.colorbar(c1)
    plt.savefig(f"res/T_{time}.png")
    plt.close()
# ========================================
# Résolution
# ========================================

T = np.zeros((mMax, nMax))
EVOL = np.zeros((mMax, nMax))

# Conditions initiales
T[:,:] = TINF(0)-10

# Boucle de résolution

time = 0 # temps physique [s]
nextLog = 0 # instant de la prochaine sauvegarde

while time<D:
    # Sauvegarde ?
    if time>=nextLog:
        export()
        nextLog += dtLog # planification de la prochaine sauvegarde
    
    # 1) Calcul du nouveau H
    for m in range(mMax):
        for n in range(nMax):

            # Calcul du flux GAUCHE
            if m==0:
                FG = 0
            else:
                ktmp = K(0.5*(T[m-1, n]+T[m, n]))
                FG = ktmp*Sx*(T[m-1, n]-T[m, n])/dx
            
            # Calcul du flux DROIT
            if m==mMax-1:
                FD = H2*Sx*(TINF(time)-T[m,n])
            else:
                ktmp = K(0.5*(T[m+1, n]+T[m, n]))
                FD = ktmp*Sx*(T[m+1,n]-T[m, n])/dx
            
            # Calcul du flux BAS
            if n==0:
                FB = 0
            else:
                ktmp = K(0.5*(T[m, n-1]+T[m, n]))
                FB = ktmp*Sy*(T[m, n-1]-T[m, n])/dy
            
            # Calcul du flux HAUT
            if n==nMax-1:
                FH = Sy*( H1*(TINF(time)-T[m,n]) + PHI)
            else:
                ktmp = K(0.5 * (T[m, n + 1] + T[m, n]))
                FH = ktmp*Sy*(T[m,n+1]-T[m, n])/dy
        
        
            # calculer le bilan
            EVOL[m, n] = dt*(FG+FD+FH+FB)/(rho*V*C(T[m,n]))
    
    # 2) préparer l'itération suivante
    T[:,:] = T[:,:] + EVOL[:,:]
    time = time + dt


