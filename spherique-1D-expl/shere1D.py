#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 15:24:32 2022

@author: sgibout
"""
import numpy as np
import matplotlib.pyplot as plt


# ========================================
# Déclaration des paramètres
# ========================================

mMax = 15       # [-] nombre de pas d'espace

R = 0.1         # [m] épaisseur du mur
D = 720        # [s] durée de la simulation

dt = 0.1        # [s] pas de temps de simulation
dtLog = 72       # [s] pas de temps de sauvegardeè2

rho = 2*240       # [kg/m3] masse volumique
C  = lambda T : 800   + 10*(T-20)       # [kg/m3] capacité calorifique
k = lambda T : 1.6 + 0.01*(T-20)         # [W/(m.K)] conductivité thermique

phi =  1000     # densité de flux imposée W/m2

T0 = 20         # [°C] température initiale



# ========================================
# préparation des géométries
# ========================================

dr = R/mMax
S = np.zeros(mMax)
V = np.zeros(mMax)

r = lambda m : (m+1-0.5)*dr

for m in range(0,mMax):
    S[m] = 4*np.pi*(r(m+0.5)**2)
    V[m] = 4*np.pi*(r(m+0.5)**3-r(m-0.5)**3)/3
    
# ========================================
# Routine d'export
# ========================================
data = np.zeros((0,3)) # variable pour le tracé de l'évolution temporelle
def export():
    global data
    x = [ r(m) for m in range(mMax)]
    plt.plot(x,T, label=f"{time:10.2f}s")
    plt.legend()
    data = np.append(data, [[time, T[0], T[mMax-1]]], axis=0)
   
# ========================================
# Résolution
# ========================================

T = np.zeros(mMax)
T[:] = T0       # Condition initiale
EVOL = np.zeros(mMax)


time = 0 # temps physique [s]

nextLog = 0 # instant de la prochaine sauvegarde

while time<D:
    # Sauvegarde ?
    if time>=nextLog:
        export()
        nextLog += dtLog # planification de la prochaine sauvegarde
    
    for m in range(mMax):
        # Calcul du flux GAUCHE
        if m==0:
            FG = 0
        else:
            k_tmp = k(0.5*(T[m-1]+T[m])) 
            FG = k_tmp*S[m-1]*(T[m-1]-T[m])/dr
        
        # Calcul du flux DROIT
        if m==mMax-1:
            FD = phi*S[m]
        else:
            k_tmp = k(0.5*(T[m+1]+T[m])) 
            FD = k_tmp*S[m]*(T[m+1]-T[m])/dr
        
        # calculer le bilan
        EVOL[m] = (FG+FD)/(rho*C(T[m])*V[m])
    
    # préparer l'itération suivante
    T[:]=T[:] + EVOL[:]
    time = time + dt
    

plt.savefig(f"graph_full.png")
plt.close()

plt.plot(data[:,0],data[:,1],label="centre")
plt.plot(data[:,0],data[:,2],label="paroi")
plt.legend()
plt.savefig("graph_evol.png")
plt.close

    