import numpy as np


def charge(filename):

    # on charge le fichier
    data = np.loadtxt(filename, delimiter=",")

    # on extrait le temps et le volume consommé sur l'intervalle
    t = data[:,0]
    v = data[:,1]

    # On vérifie que le pas est le bon (cf. CAS1)
    tmp = t[1:]-t[0:-1]
    dt = np.mean(tmp)
    if np.std(tmp)>dt/1000: # par exemple:
        raise Exception("Pas de temps non régulier dans le profil")

    return dt, len(v), v/dt

