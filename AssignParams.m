%%  Set input parameters here, then save and run RunSimulation.m

function [planet,t,orbInit,orbFin,numHT,numLT,sc,prop] = AssignParams(thrust_level)

%% SPECIFY PLANETARY PROPERTIES
Re = 6378.15; %Radius of Earth in km
omega_E = 7.292e-5; %rotation rate of the earth in rad/s
day_sd = 86164.1; %number of seconds in a sidereal day
mu = 398600; %Grav. param. for Earth in km^3/s^2
g0 = 9.81; % Gravitational acceleration in m/s^2;
J2 = 1.08263e-3; 

%% SPECIFY SIMULATION TIME VECTOR
% t = [1:5*60:2*3600*24]; %Time vector in s - take position every 5 mins
numDays = 50; %number of days you're looking to simulate
numPoints = 10000; %number of data points within the simulation time span
        
        t = linspace(0,numDays*day_sd, numPoints);%linearly spaced tvec points in s (better for analysis w/o plotting, no periodic behavior). 

%% SPECIFY INITIAL ORBIT ELEMENTS
%ROAMS GOM orbit parameters (from 16.851 Fall 2021 Final Pres.)
initial_alt = 497.5; %initial altitude in km
e = 0;
incl = 51.6; % inclination in deg
RAAN = 60; % in deg. 
ArgPer = 0; % in deg. initially say 0 for placeholder
anomaly = 0; % in deg. initially say 0 for placeholder

        a = initial_alt+Re;
        h = sqrt(a*mu);
        T = 2*pi/(sqrt(mu))*a^(3/2); %Orbital period in s

%% SPECIFY FINAL DESIRED ORBIT ELEMENTS
RAAN_f = 59; %in deg. 
    final_alt = initial_alt; %in km, should be same as initial_alt for ROM->ROM transfers
alt_tol = 0.2; %Tolerance for final altitude [km]. Default is 0.2 km. %NOTE: May not see this work well for large timesteps

GOM_alt_diff = 50; % ROM-GOM altitude difference in km (if simulating GOM transfers). We can go +diff or -diff to move in different RAAN directions.
        
        delta_RAAN = RAAN_f-RAAN; %difference in deg. NOTE: Revisit for changeover back to 0 (prevent differencing issue). 
        r_0 = a; %initial orbital radius
        r_f = final_alt + Re;

%% SPECIFY CONSTRAINTS ON MANEUVERS
% t_f = 10;%Time in DAYS
% t_f = t_f*day_sd; %convert to seconds

numHT = 2; %number of desired ROM->ROM maneuvers with high thrust (will be used to determine elliptical transfer properties)
numLT = 5; %number of desired ROM->ROM maneuvers with low thrust 

%% SPECIFY SPACECRAFT PARAMETERS
m_0 = 12; %spacecraft wet mass in kg

%% SPECIFY PROPULSION SYSTEM PARAMETERS
if thrust_level == 1 %HIGH THRUST
    thrust = 1; %nominal thrust for Aeroject GR1 in N 
    I_sp = 220; % (May not be the correct number for Aerojet GR1)
    prop_mass = 0.5; % in kg for 1U Aerojet GR1
    prop_power = 20; %Power in W (assumed 20 - not necessarily accurate for GR1)
    
elseif thrust_level == 2 %LOW THRUST
    thrust = 330e-6; %nominal thrust for Enpulsion Nano in N (https://www.enpulsion.com/wp-content/uploads/ENP2018-001.G-ENPULSION-NANO-Product-Overview.pdf)
    I_sp = 6000; %specific impulse in s for impulsion Nano (high estimate)
    prop_power = 40;% Power of electric propulsion system in W
    prop_mass = 0.220; % initial mass of propellant in kg
end

    accel = (thrust/m_0)/1000; %in km/s
    m_dot = thrust/(I_sp*g0); % mass flow rate in kg/s;
    t_thrust = prop_mass/m_dot; % total time in s that thruster can burn before running out of propellant
    dv_avail = t_thrust*accel; %dv available in km/s
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DO NOT EDIT BELOW
%% SPECIFY PLANETARY PROPERTIES
planet.Re = Re; 
planet.omega_E = omega_E;
planet.day_sd = day_sd;
planet.mu = mu;
planet.g0 = g0;
planet.J2 = J2;

%% SPECIFY INITIAL ORBIT ELEMENTS
%ROAMS GOM orbit parameters (from 16.851 Fall 2021 Final Pres.)
orbInit.initial_alt = initial_alt;
orbInit.e = e;
orbInit.incl = incl;
orbInit.RAAN = RAAN;
orbInit.ArgPer = ArgPer;
orbInit.anomaly = anomaly;

        orbInit.a = a;
        orbInit.h = h;
        orbInit.T = T;
        
%% SPECIFY FINAL DESIRED ORBIT ELEMENTS
orbFin.RAAN_f = RAAN_f;
    orbFin.final_alt = final_alt; 
orbFin.alt_tol = alt_tol;
orbFin.GOM_alt_diff = GOM_alt_diff;
        
        orbFin.delta_RAAN = delta_RAAN; 
        orbFin.r_0 = r_0;
        orbFin.r_f = r_f;

%% SPECIFY SPACECRAFT PARAMETERS
sc.m_0 = m_0;

%% SPECIFY PROPULSION SYSTEM PARAMETERS
if thrust_level == 1 %HIGH THRUST
    prop.thrust = 1; %nominal thrust for Aeroject GR1 in N 
    prop.I_sp = 235; % (May not be the correct number for Aerojet GR1)
    prop.prop_mass = 0.5; % in kg for 1U Aerojet GR1
    prop.prop_power = 20; %Power in W (assumed 20 - not necessarily accurate for GR1)
    
elseif thrust_level == 2 %LOW THRUST
    prop.thrust = 330e-6; %nominal thrust for Enpulsion Nano in N (https://www.enpulsion.com/wp-content/uploads/ENP2018-001.G-ENPULSION-NANO-Product-Overview.pdf)
    prop.I_sp = 3500; %specific impulse in s for impulsion Nano (high estimate)
    prop.prop_power = 40;% Power of electric propulsion system in W
    prop.prop_mass = 0.220; % initial mass of propellant in kg
end

    sc.accel = accel;
    sc.m_dot = m_dot;
    sc.t_thrust = thrust;
    sc.dv_avail = dv_avail;

end
