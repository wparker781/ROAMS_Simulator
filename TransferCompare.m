%% Function Description
% 
% Run this function from Run_TransferCompare.m. In that script, you can
% iterate between transfer types and the desired change in RAAN to help
% determine which transfer type is best for a given new target location. 
%
% Developed by William Parker - March 2021

function [t_trans, dv, num_trans] = TransferCompare(trans_type, thrust_level, delta_RAAN)

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
RAAN = 60; % in deg. initially say 0 for placeholder
ArgPer = 0; % in deg. initially say 0 for placeholder
anomaly = 0; % in deg. initially say 0 for placeholder

    a = initial_alt+Re;
    h = sqrt(a*mu);
    T = 2*pi/(sqrt(mu))*a^(3/2); %Orbital period in s

%% SPECIFY FINAL DESIRED ORBIT ELEMENTS
final_alt = 497.5; %in km, should be same as initial_alt for ROM->ROM transfers
    % RAAN_f = 70; %in deg. 
    % delta_RAAN = RAAN_f-RAAN; %difference in deg. NOTE: Revisit for changeover back to 0 (prevent differencing issue). 
alt_tol = 0.2; %Tolerance for final altitude [km]. Default is 0.2 km. %NOTE: May not see this work well for large timesteps

GOM_alt_diff = 50; % ROM-GOM altitude difference in km (if simulating GOM transfers). We can go +diff or -diff to move in different RAAN directions.

    r_0 = a; %initial orbital radius
    r_f = final_alt + Re;

%% SPECIFY CONSTRAINTS ON MANEUVERS
% t_f = 10;%Time in DAYS
% t_f = t_f*day_sd; %convert to seconds

numHT = 2; %number of desired ROM->ROM maneuvers with high thrust
numLT = 5; %number of desired ROM->ROM maneuvers with low thrust

%% SPECIFY SPACECRAFT PARAMETERS
m_0 = 12; %spacecraft wet mass in kg

%% SPECIFY PROPULSION SYSTEM PARAMETERS
if thrust_level == 1
    %HIGH THRUST
    thrust = 1; %nominal thrust for Aeroject GR1 in N 
    I_sp = 220; % (May not be the correct number for Aerojet GR1)
    prop_mass = 0.5; % in kg for 1U Aerojet GR1
    prop_power = 20; %Power in W (assumed 20 - not necessarily accurate for GR1)
elseif thrust_level == 2
    %LOW THRUST
    thrust = 330e-6; %nominal thrust for Enpulsion Nano in N (https://www.enpulsion.com/wp-content/uploads/ENP2018-001.G-ENPULSION-NANO-Product-Overview.pdf)
    I_sp = 6000; %specific impulse in s for impulsion Nano (high estimate)
    prop_power = 40;% Power of electric propulsion system in W
    prop_mass = 0.220; % initial mass of propellant in kg
end

accel = thrust/m_0;
m_dot = thrust/(I_sp*g0); % mass flow rate in kg/s;
t_thrust = prop_mass/m_dot; % total time in s that thruster can burn before running out of propellant
dv_avail = t_thrust*accel/1000; %dv available in km/s

% Note: For the current implementation, we assume impulsive maneuvers and
% constant mass, so specifics for the high thrust system are only used to
% determine maximum delta-v. 

%% PROPAGATE SCENARIO WITH J2 CORRECTION
% At the moment, we're only able to propagate low thrust transfers (high
% thrust are done from analytical approximations, not numerical
% simulation).


%if high thrust... analytically determine transfer time, dv requirement,
%number of transfers possible

if trans_type == 1 %GOM transfer
    [t_trans, dv, num_trans] = GOMTrans(GOM_alt_diff, initial_alt, final_alt, dv_avail, e, incl, delta_RAAN, Re, mu, J2);
    fprintf('<strong>Circular GOM transfer completed (using high thrust)!</strong>\n')

elseif trans_type == 2 %Elliptical transfer
    %assign dv we're willing to sacrifice per maneuver (>dv, <t)
    dv_per_trans = dv_avail/numHT; %divide total dv by number of transfers we'd like
    [t_trans, dv, num_trans] = ellipticalTrans(initial_alt, dv_avail,e,incl, delta_RAAN, Re, mu, J2, dv_per_trans);
    fprintf('<strong>Elliptical transfer completed (using high thrust)!</strong>\n')

elseif trans_type == 3 %low thrust spiral transfer
    RAAN_f = RAAN+delta_RAAN;
    [posN,RAAN, RAAN_dot, RAAN_rel, ArgPer,v,alt,thrust, t, t_trans, dv, num_trans] = spiralTrans(t,mu,a,e,h,incl,RAAN,ArgPer,anomaly, thrust, m_0,m_dot,prop_power, r_f, alt_tol, RAAN_f, Re, dv_avail);
    fprintf('<strong>Time-optimal spiral transfer completed (using low thrust)!</strong>\n')

elseif trans_type == 4 %Direct Plane Change
    [t_trans, dv, num_trans] = PlaneChangeTrans(initial_alt, dv_avail,e,incl, delta_RAAN, Re, mu, J2);
    fprintf('<strong>Direct plane change transfer completed (using high thrust)!</strong>\n')
end

fprintf('<strong>Transfer time:</strong> %.2f seconds (or %.3f hours or %.3f days) \n', t_trans, t_trans/3600, t_trans/(day_sd));
fprintf('<strong>Delta-v required:</strong> %.4f km/s or %.2f m/s \n', dv, dv*1000);
fprintf('<strong>Number of transfers possible with this prop system and transfer type:</strong> %.3f \n', num_trans);



end


