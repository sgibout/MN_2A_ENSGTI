#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 15:24:32 2022

@author: sgibout
"""
import numpy as np
import matplotlib.pyplot as plt
import sys


# ========================================
# Déclaration des paramètres
# ========================================

mMax = 30       # [-] nombre de pas d'espace selon X
nMax = 30       # [-] nombre de pas d'espace selon Y

LX = 0.1         # [m] longueur mur selon X
LY = 0.12         # [m] longueur mur selon Y
D = 600        # [s] durée de la simulation

dt = 15        # [s] pas de temps de simulation
dtLog = 10       # [s] pas de temps de sauvegardeè2

rho = 1000       # [kg/m3] masse volumique
CL= 4000   # [J/kg/K] capacité calorifique état liquide
CS= 2000   # [J/kg/K] capacité calorifique état solide
TF = 30     # [°C] température de fusion
LF = 100E3  # [J/kg] chaleur latente de fusion
k = 1         # [W/(m.K)] conductivité thermique

phi =  10000     # densité de flux imposée W/m2

T0 = 20         # [°C] température initiale
TINF = 50       # [°C] température fluide extérieur
HCONV = 100         # [W/m2/K] Coefficient d'échange global


# ========================================
# préparation des géométries
# ========================================

dx = LX/mMax
dy = LY/nMax
Sx = 1*dy
Sy = 1*dx
V = 1*dx*dy


Hx = 1/(1/HCONV + 0.5*dx/k)
Hy = 1/(1/HCONV + 0.5*dy/k)


# ========================================
# Routine d'export
# ========================================
def export():
#    plt.figure()
    c1=plt.contourf(T,     cmap ="jet")
    plt.colorbar(c1)
    plt.savefig(f"res/T_{time}.png")
    plt.close()
#    plt.figure()
    c1=plt.contourf(Y,     cmap ="jet")
    plt.colorbar(c1)
    plt.savefig(f"res/Y_{time}.png")
    plt.close()    
# ========================================
# Résolution
# ========================================

H = np.zeros((mMax, nMax))
T = np.zeros((mMax, nMax))
Y = np.zeros((mMax, nMax))

# Conditions initiales
T[:,:]= T0
if T0<TF:
    # SOLIDE
    H[:,:] = (T0-TF)*CS
    Y[:,:] = 0
elif T0>TF:
    # Liquide
    H[:,:] = LF+(T0-TF)*CL
    Y[:,:] = 1
else:
    print("*** ERREUR : impossible d'initialiser à T0=TF")
    sys.exist()
    


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
                FG = k*Sx*(T[m-1, n]-T[m, n])/dx
            
            # Calcul du flux DROIT
            if m==mMax-1:
                FD = Hx*Sx*(TINF-T[m,n])
            else:
                FD = k*Sx*(T[m+1,n]-T[m, n])/dx
            
            # Calcul du flux BAS
            if n==0:
                FB = phi*Sy
            else:
                FB = k*Sy*(T[m, n-1]-T[m, n])/dy
            
            # Calcul du flux HAUT
            if n==nMax-1:
                FH = Hy*Sy*(TINF-T[m,n])
            else:
                FH = k*Sy*(T[m,n+1]-T[m, n])/dy
        
        
            # calculer le bilan
            H[m, n] += dt*(FG+FD+FH+FB)/(rho*V)
    
    # 2) Evolution vers le nouvel état (mise à jour de T et Y)
    for m in range(mMax):
        for n in range(nMax):
            if H[m,n]<0:
                # Solide
                T[m,n] = TF+H[m,n]/CS
                Y[m,n] = 0
            elif H[m,n]> LF:
                # Liquide
                T[m,n] = TF+(H[m,n]-LF)/CL
                Y[m,n] = 1
            else:
                # Equilibre
                T[m,n] = TF
                Y[m,n] = H[m,n]/LF
                
    
                
    
    
    # 3) préparer l'itération suivante
    time = time + dt
    
