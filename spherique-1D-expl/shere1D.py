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

R = 0.1         # [m] rayon de la sphère
D = 720        # [s] durée de la simulation

dt = 0.1        # [s] pas de temps de simulation
dtLog = 1       # [s] pas de temps de sauvegardeè2

rho = 2*240       # [kg/m3] masse volumique
C  = 800           # [J/kg/K] capacité calorifique
k = 0.5 # [W/(m.K)] conductivité thermique

epsilon =  0.7     # émissivité de la paroi [-]
Trad = 35       # Température radiative l'ambiance externe [°C]

H = 50         # Coefficient d'échange convectif à la paroi [W/m2/K]
Tinf = 25       # Température air externe [°C]


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
data = np.zeros((0,4)) # variable pour le tracé de l'évolution temporelle
def export():
    global data
    # On mémorise les données qu'on soutaite manipuler (ici centre et dernier noeud et paroi)
    data = np.append(data, [[time, T[0], T[mMax-1], TP]], axis=0)

def calcTP(TM,TP):
    # On utilie une méthode Newton classique
    f = lambda x : 2*k*(TM-x)/dr + H*(Tinf-x) + epsilon*5.67E-8*((Trad+273.15)**4-(x+273.15)**4)
    df = lambda  x : -2*k/dr - H - 4*epsilon*5.67E-8*(x+273.15)**3

    delta = np.inf
    x = TP
    while abs(delta)>1e-3:
        delta = f(x)/df(x)
        x = x - delta

    return x
# ========================================
# Résolution
# ========================================

T = np.zeros(mMax)
TP = T[:] = T0       # Condition initiale
EVOL = np.zeros(mMax)

Uinit = np.sum(V*T)*rho*C
sumFlux = 0
time = 0 # temps physique [s]

nextLog = 0 # instant de la prochaine sauvegarde

while time<D:
    # Sauvegarde ?
    if time>=nextLog:
        export()
        nextLog += dtLog # planification de la prochaine sauvegarde

    # Calcul de la température de paroi TP
    TP = calcTP(T[-1], TP) # On calcule le nouveau TP (Newton) en partant de TP comme estimation initiale

    for m in range(mMax):
        # Calcul du flux GAUCHE
        if m==0:
            FG = 0
        else:
            FG = k*S[m-1]*(T[m-1]-T[m])/dr
        
        # Calcul du flux DROIT
        if m==mMax-1:
            FD = S[m]* (H*(Tinf - TP) + epsilon*5.67E-8*((Trad+273.15)**4-(TP+273.15)**4) )
            sumFlux += FD
        else:
            FD = k*S[m]*(T[m+1]-T[m])/dr
        
        # calculer le bilan
        EVOL[m] = (FG+FD)/(rho*C*V[m])
    
    # préparer l'itération suivante
    T[:]=T[:] + EVOL[:]
    time = time + dt
    
Ufinal = np.sum(V*T)*rho*C

# On trace la figure à partir des données mémorisées au cours de la résolution
plt.plot(data[:,0],data[:,1],label="Noeud central")
plt.plot(data[:,0],data[:,2],label="Noeud périphérie")
plt.plot(data[:,0],data[:,3],label="Paroi")
plt.xlabel('t [s]')
plt.ylabel('T(t) [°C]')
plt.legend()
plt.savefig("graph_evol.png")
plt.close()

print("Vérification régime permanent:")
print(f"densité de flux convectif : {H*(T[-1]-Tinf)}")
print(f"densité de flux radiatif : {epsilon*5.67E-8*((T[-1]+273.15)**4-(Trad+273.15)**4)}")

print("Vérification conservation:")
print(f"Ufinal-Uinit = {Ufinal-Uinit}")
print(f"SomFlux      = {sumFlux}")
