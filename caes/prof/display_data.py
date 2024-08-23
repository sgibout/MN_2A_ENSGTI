import numpy as np
from matplotlib import pyplot as plt


# Données brutes à la seconde sur 365 jours
GHI = np.loadtxt('GHI_année_1 mn.txt')
conso = 1000*np.loadtxt('semaine_conso.txt') # en kW

t= np.arange(len(GHI))/60/24 # t est donc en seconde

# Données à la seconde sur 7 jours
conso_r = np.tile(conso,53)[:len(GHI)]
plt.plot(t,conso_r,'r')
#plt.plot(t,GHI)
#plt.show()

data = np.column_stack((t,GHI,conso_r))
plt.plot(data[:,0],data[:,1])
plt.show()

np.save('../data.npy', data)


#conso_r = np.tile(conso,53)

#t_CONSO = np.arange(len(conso))/60/24 # t est donc en seconde
#t_CONSO_r = np.tile(t_CONSO,53)

#plt.plot(t_CONSO_r,conso_r)
#plt.plot(t_GHI,GHI,"r")
#plt.show()

