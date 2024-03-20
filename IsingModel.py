import numpy as np
import matplotlib.pyplot as plt
import random

#=======================================================
#here you describe the parameters of the system
Nx = 20 #number of spins in x direction
Ny = 20 #number of spins in y direction
J = 1 #interaction energy
kB = 1 #Boltzman constant
Tc = 2 * J / (kB * np.log(np.sqrt(2) + 1)) #critical temperature
nSteps = 1000000 #number of steps for Metropolis algorithm
#=======================================================

#=======================================================
#all necessary functions
def Spin_Configuration(Nx, Ny): #defines spin configuration
  return np.random.choice([-1, 1], size=(Ny, Nx))

def Plot(S): #plot of spin configuration
  x = np.arange(0, Nx, 1)
  y = np.arange(0, Ny, 1)

  fig, ax = plt.subplots(1, 1)

  c = ax.pcolormesh(x, y, S)
  fig.colorbar(c, ax=ax)
  plt.show()

def Energy(S): #calculates energy of spin configuration
  H = 0

  for i in range(Ny-1):
    for j in range(Nx-1):
      H = H - J * (S[i][j] * S[i][j+1] + S[i][j] * S[i+1][j])

  return H

def Magnetization(S): #calculates magnetization of spin configuration
  M = np.sum(S)
  return M

def Energy_S(i,j): #calculates interaction energy for element ij
  if i + 1 == Ny:
    i_next = 0
  else:
    i_next = i + 1

  if i - 1 < 0:
    i_prev = Ny - 1
  else:
    i_prev = i - 1

  if j + 1 == Nx:
    j_next = 0
  else:
    j_next = j + 1

  if j - 1 < 0:
    j_prev = Nx - 1
  else:
    j_prev = j - 1

  Hij = - J * (S[i][j] * S[i_next][j] + S[i][j] * S[i_prev][j] +
              S[i][j] * S[i][j_next] + S[i][j] * S[i][j_prev])
  return Hij

def metropolis(nSteps,T): #metropolis algorithm
  Mobs = []
  Eobs = []

  Eobs.append(Energy(S))
  Mobs.append(Magnetization(S))

  for k in range(1,nSteps):
    i = random.randrange(0,Nx,1)
    j = random.randrange(0,Ny,1)
    E0 = Energy_S(i,j)
    S[i][j] = - S[i][j]
    E1 = Energy_S(i,j)
    dE = E1 - E0

    if dE < 0:
      Eobs.append(Eobs[k-1] + dE)
      if S[i][j] == -1:
        Mobs.append(Mobs[k-1] - 2)
      else:
        Mobs.append(Mobs[k-1] + 2)
    elif (dE > 0 and np.random.uniform(0,1) > np.exp(-dE / T)):
      S[i][j] = - S[i][j]
      Eobs.append(Eobs[k-1])
      Mobs.append(Mobs[k-1])
    else:
      Eobs.append(Eobs[k-1] + dE)
      if S[i][j] == -1:
        Mobs.append(Mobs[k-1] - 2)
      else:
        Mobs.append(Mobs[k-1] + 2)

  return Eobs, Mobs

def Observable_Magnetisation(Mobs,T): #avarage magnetization
  M = 0
  
  for i in range(int(len(Mobs) * 0.5)):
    M = M + Mobs[i + int(len(Mobs) * 0.5)] / int(len(Mobs) * 0.5)

  return abs(M)

def Heat_capacity(Eobs,T): #heat capacity
  E = 0
  E2 = 0

  for i in range(int(len(Eobs) * 0.5)):
    E = E + Eobs[i + int(len(Eobs) * 0.5)] / int(len(Eobs) * 0.5)
    E2 = E2 + (Eobs[i + int(len(Eobs) * 0.5)] ** 2) / int(len(Eobs) * 0.5)

  C = kB * (E2 - E ** 2) / ((kB * T) ** 2)

  return C

def Observable_Energy(Eobs, T): #energy
  E = 0
  
  for i in range(int(len(Eobs) * 0.5)):
    E = E + Eobs[i + int(len(Eobs) * 0.5)] / int(len(Eobs) * 0.5)

  return E

#=======================================================

#=======================================================
#here you can calculate what do you want


#lets show spin configuration for different temperatures
T = [0.5 * Tc, 1 * Tc, 2 * Tc]

for t in T:
  print('Stage', np.where(T == t)[0][0],', T=',t/Tc,'Tc :started')
  S = Spin_Configuration(Nx, Ny)
  Plot(S)
  [Eobs,Mobs] = metropolis(nSteps,t)
  Plot(S)
  print('Stage', np.where(T == t)[0][0],', T=',t/Tc,'Tc :finished')

#lets calculate temperature dependencies for magnetization and heat capacity
T = np.linspace(0.7 * Tc, 3 * Tc, 50)
M = []
C = []
E = []

for t in T:
  print('Stage', np.where(T == t)[0][0],', T=',t/Tc,'Tc :started')
  m = 0
  c = 0
  e = 0
  for k in range(50):
    S = Spin_Configuration(Nx, Ny)
    [Eobs,Mobs] = metropolis(nSteps,t)
    m = m + Observable_Magnetisation(Mobs,t)
    c = c + Heat_capacity(Eobs,t)
    e = e + Observable_Energy(Eobs,t)

  M.append(m / 50)
  C.append(c / 50)
  E.append(e / 50)
  print('Stage', np.where(T == t)[0][0],', T=',t/Tc,'Tc :finished')

plt.plot(T/Tc, M)
plt.xlabel('T/Tc')
plt.ylabel('M')
plt.show()

plt.plot(T/Tc, C)
plt.xlabel('T/Tc')
plt.ylabel('C')
plt.show()

plt.plot(T/Tc, E)
plt.xlabel('T/Tc')
plt.ylabel('E')
plt.show()
