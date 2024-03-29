%Script for testing PropagateOrbit_J2 for ROAMS
% William Parker - March 2021
clc; clear all; close all; 

%% SPECIFY EARTH PROPERTIES
Re = 6378.15; %Radius of Earth in km
omega_E = 7.292e-5; %rotation rate of the earth in rad/s
day_sd = 86164.1; %number of seconds in a sidereal day
mu = 398600; %Grav. param. for Earth in km^3/s^2

%% SPECIFY SIMULATION TIME VECTOR
% t = [1:5*60:2*3600*24]; %Time vector in s - take position every 5 mins
t = linspace(0,10*day_sd, 100);%linearly spaced tvec points in s (better for analysis w/o plotting, no periodic behavior). 

%% SPECIFY INITIAL ORBIT PARAMETERS
%ROAMS GOM orbit parameters (from 16.851 Fall 2021 Final Pres.)
initial_alt = 547.5; %initial altitude in km
a = initial_alt+Re;
e= 0;
h = sqrt(a*mu);
incl = 51.6; % inclination in deg
RAAN = 0; %initially say 0 for placeholder
ArgPer = 0; %initially say 0 for placeholder
anomaly = 0; %initially say 0 for placeholder

T = 2*pi/(sqrt(mu))*a^(3/2); %Orbital period in s

%% SPECIFY FINAL DESIRED ORBIT PARAMETERS
final_alt = 497.5;
alt_tol = 0.2; %Tolerance for final altitude [km]. Default is 0.2 km. %NOTE: May not see this work well for large timesteps

r_0 = a; %initial orbital radius
r_f = final_alt + Re;

%% SPECIFY SPACECRAFT PARAMETERS
m_0 = 8; %spacecraft wet mass in kg

%% SPECIFY PROPULSION SYSTEM PARAMETERS
max_thrust = 330e-6; %nominal thrust for Enpulsion Nano (https://www.enpulsion.com/wp-content/uploads/ENP2018-001.G-ENPULSION-NANO-Product-Overview.pdf)
I_sp = 2000; %specific impulse in s
prop_power = 40;% Power of electric propulsion system in W
m_dot = max_thrust/(9.81*I_sp); % mass flow rate from prop %%!!! Modify for actual thrust in propellant consumtion


[posN,RAAN, RAAN_dot ,v,alt,thrust, t] = PropagateOrbit_J2(t,mu,a,e,h,incl,RAAN,ArgPer,anomaly, max_thrust, m_0,m_dot,prop_power, r_f, alt_tol, Re);


%% Plot Results

% Plot sphere with 3D scatterplot of orbital locations (Doesn't make much
% sense unless the timestep is <5m). 
    % figure()
    % %make and plot sphere for earth
    % [sx,sy,sz] = sphere();
    % r = 6378;
    % surf( r*sx, r*sy, r*sz )
    % hold on
    % axis equal
    % %plot propagated orbit at times t
    % plot3(posN(1,:),posN(2,:),posN(3,:),'r.');
    % xlabel('x [km]');
    % ylabel('y [km]');
    % zlabel('z [km]');
    % hold off

% Plot ratio of earth's rotation rate to RAAN precession rate (when it
% reaches an integer value, you've reached an RGT orbit). 
    figure()
    plot(alt, rad2deg(omega_E)./rad2deg(RAAN_dot));
    xlabel('Altitude [km]')
    ylabel('Omega_E/RAAN dot');
    %ideally, we want the rotation rate of the earth to be an integer multiple of
    %the rate of precession of the orbit
    hold off

%Plot RAAN of orbit over time
    figure()
    plot(t, RAAN)
    xlabel('Time [s]')
    ylabel('Orbital Longitudinal Precession(RAAN)[deg]')

%Plot altitude over time
    figure()
    plot(t, alt)
    xlabel('Time [s]')
    ylabel('Altitude [km]');

%Plot Thrust over time
    figure()
    plot(t, thrust)
    xlabel('Time [s]')
    ylabel('Commanded Thrust [N]')



