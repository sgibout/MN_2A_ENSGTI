#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: sgibout
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# ========================================
# Déclaration des paramètres
# ========================================

mMax = 30       # [-] nombre de pas d'espace selon X
nMax = 30       # [-] nombre de pas d'espace selon Y

LX = 0.05         # [m] longueur mur selon X
LY = 0.05         # [m] longueur mur selon Y
D = 3600        # [s] durée de la simulation

dt = 0.1       # [s] pas de temps de simulation
dtLog = 60       # [s] pas de temps de sauvegardeè2

rho = 1000       # [kg/m3] masse volumique

CL= 4000    # [J/kg/K] capacité calorifique état liquide
CS= 2000    # [J/kg/K] capacité calorifique état solide

KL = 2      # [W/m/K] conductivité thermique état liquide
KS = 3      # [W/m/K] conductivité thermique état solide

TF = 35     # [°C] température de fusion
LF = 100E3  # [J/kg] chaleur latente de fusion

TINF = lambda t : 30 +  10*np.sin(2*np.pi*t/3600) # [°C] température externe

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
XX, YY = np.meshgrid(
    [ (m+0.5)*dx for m in range(0,mMax) ],
    [ (n+0.5)*dy for n in range(0,nMax) ]
)

def export():
    Tlevels = np.linspace(20, 55, num=10)
    fig,(axT,axY) = plt.subplots(1,2,figsize=(18,6))
    c1 = axT.pcolormesh(XX,YY,T,     cmap ="jet")
    c2 = axY.pcolormesh(XX,YY,Y,     cmap ="jet")
    axT.x_label="x [m]"
    axT.y_label="y [m]"
    axY.x_label="x [m]"

    cbar1 = fig.colorbar(c1, ax=axT)
    cbar2 = fig.colorbar(c2, ax=axY)

    plt.suptitle(f"t={time:.0f}s")
    plt.savefig(f"res/T_{time}.png")
    plt.close()
# ========================================
# Résolution
# ========================================

U = np.zeros((mMax, nMax))
T = np.zeros((mMax, nMax))
Y = np.zeros((mMax, nMax))


# Conditions initiales
# Conditions initiales
T[:,:] = T0 = TINF(0)
if T0<TF:
    # SOLIDE
    U[:,:] = (T0-TF)*CS
    Y[:,:] = 0
elif T0>TF:
    # Liquide
    U[:,:] = LF+(T0-TF)*CL
    Y[:,:] = 1
else:
    print("*** ERREUR : impossible d'initialiser à T0=TF")
    sys.exit(-1)
# Boucle de résolution

time = 0 # temps physique [s]
nextLog = 0 # instant de la prochaine sauvegarde

# fonction utilitaire pour le calcul de K à l'interface entre deux nœuds qui
# possède les fractions liquides Y1 et Y2
calcK = lambda Y1,Y2 : KS + 0.5*(Y1+Y1)*(KL-KS)

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
                FG = calcK(Y[m-1,n],Y[m,n])*Sx*(T[m-1, n]-T[m, n])/dx
            
            # Calcul du flux DROIT
            if m==mMax-1:
                FD = H2*Sx*(TINF(time)-T[m,n])
            else:
                FD = calcK(Y[m+1,n],Y[m,n])*Sx*(T[m+1,n]-T[m, n])/dx
            
            # Calcul du flux BAS
            if n==0:
                FB = 0
            else:
                FB = calcK(Y[m,n-1],Y[m,n])*Sy*(T[m, n-1]-T[m, n])/dy
            
            # Calcul du flux HAUT
            if n==nMax-1:
                FH = Sy*( H1*(TINF(time)-T[m,n]) + PHI)
            else:
                FH = calcK(Y[m,n+1],Y[m,n])*Sy*(T[m,n+1]-T[m, n])/dy
        
        
            # calculer le bilan
            U[m, n] += dt*(FG+FD+FH+FB)/(rho*V)
    
    # 2) Evolution vers le nouvel état (mise à jour de T et Y)
    # On pourrait (et ce serait mieux) faire une fonction...
    for m in range(mMax):
        for n in range(nMax):
            if U[m, n] < 0:
                # Solide
                T[m, n] = TF + U[m, n] / CS
                Y[m, n] = 0
            elif U[m, n] > LF:
                # Liquide
                T[m, n] = TF + (U[m, n] - LF) / CL
                Y[m, n] = 1
            else:
                # Equilibre
                T[m, n] = TF
                Y[m, n] = U[m, n] / LF

    # 3) on passe au pas de temps suivant!!
    time = time + dt


