# ROAMS_Simulator

To run the ROAMS simulator, see runSimulation.m. 
In this script, input the planetary, orbit, spacecraft, and propulsion information. 
Then, that information is fed to PropagateOrbit_J2.m.
Thrust control is done using AltController_simple.m. 

For best results, use a small simulation timestep (~5-10min). For fast results, use linspace over simulation duration, and use 100 pts. 
