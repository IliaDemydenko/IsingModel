This code simulates 2d Ising model of feromagnetic without external magnetic field. It is an university project.

Functions:
- Spin_Configuration(Nx, Ny) - generates random spin configuration with sizes Nx, Ny;
- Plot(S) - plots spin configuration S;
- Energy(S) - calculates energy of configuration S;
- Magnetization(S) - calculates magnetization of configuration S;
- Energy_S(i,j) - calculates interaction energy for element ij;
- metropolis(nSteps,T) - realize implements Metropolis-Hastings algorithm for nSteps and temperature T;
- Probability(Eobs,T) - calculates probability of each configuration based on Gibbs distribution;
- Observable_Magnetisation(Eobs,Mobs,T) - calculates observable magnetization;
- Heat_capacity(Eobs,T) - calculates observable heat capacity.

NOTE:
1) All parameters of sistem you can set at the beginning.
2) For calculations of the observable parameters we use only the last half of massives Eobs and Mobs, due to the thermalization.
3) At the beginning you can specify units by setting k and J. In the basic version they set as 1.
4) In the last section you can calculate what do you want. In basic version it visualize some spin configurations and calculates magnetization and heat capacity dependences on temperature.
