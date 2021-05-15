# ROAMS_Simulator

WELCOME TO THE ROAMS SIMULATION TOOL (RST AKA "RUSTY")

Use this tool to run ROAMS maneuver simulations. First, input planetary properties, 
simulation timing, current & desired orbit elements, spacecraft
properties, and propulsion system parameters in AssignParams.m. The, run
runSimulation.m to generate information for various methods of orbital transfer, 
and produce performance tracking metrics in the process for analysis. To
compare multiple transfer methods over a variety of transfer distances,
use Run_TransferCompare.m.

NOTE: The only files you should need to edit/run are AssignParams.m, RunSimulation.m, and Run_TransferCompare.m!

DIRECTIONS FOR USE: 
1. Edit the input parameters in AssignParams.m for desired spacecraft/orbit
2. Run RunSimulation.m.
3. Answer any questions asked in input dialog boxes or in the command
window (HT/LT = High/Low Thrust) -- These may not pop up automatically,
so if the simulation is paused for an extended period of time, look in
your matlab tabs for a selection menu window. 
4. Wait for orbit propagator to work through the entire simulation (it
will continue for the simulation duration you assigned, and not stop
automatically when the transfer has been completed). 
5. View output information for selected transfer type (should be fast to compute for
most transfer types)
6. NOTE: If you're simulating a spiral transfer and get an error that
says "Increase the duration of your simulation window," increase numDays in AssignParams.m.
If your spiral transfer plots don't appear to perform as expected, try 
increasing numPoints. The spiral transfer simulation is especially
susceptible to quantization error.
7. If you would like to compare different transfer methods over multiple
transfer scenarios (distance between targets is analagous to delta_RAAN),
edit the parameters in AssignParams.m and then run
Run_TransferCompare.m.


ASSUMPTIONS: 
1. Only circular orbits in spiral transfers
2. No inclination change in transfer, only altitude changes (except for direct plane change maneuver). 
3. Constant mass 
4. Perfect burns - exactly prograde/retrograde at nominal thrust
5. No environmental disturbance forces (drag, SRP, etc.)
6. Impulsive (instantaneous) burns for high thrust transfers.
7. No RAAN drift during high-thrust transfers to/from GOM (only while in GOM). 

CURRENT LIMITATIONS:
1. Orbit propagator is only providing data for low thrust transfers. High
thrust transfers are simulated using analytical approximations, rather
than numerical iteration through time (for impulsive maneuvers, the
maneuver itself is assumed to occur over a very small timescale)
2. Assuming purely impulsive high thrust transfers

For best results, use a small simulation timestep (~5-10min). For fast results, use linspace over simulation duration, and use 100 pts. 
