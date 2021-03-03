%Script for testing PropagateOrbit_J2 for ROAMS
% William Parker - March 2021
clc; clear all; close all; 

Re = 6378.15; %km
omega_E = 7.292e-5; %rotation rate of the earth in rad/s
day_sd = 86164.1; %number of seconds in a sidereal day

%ROAMS GOM orbit parameters (from 16.851 Fall 2021 Final Pres.)
mu = 398600; 
a = 547.5+Re;
e= 0;
h = sqrt(a*mu);
incl = 51.6; %deg
RAAN = 0; %initially say 0 for placeholder
ArgPer = 0; %initially say 0 for placeholder

anomaly = 0; %initially say 0 for placeholder

T = 2*pi/(sqrt(mu))*a^(3/2);
% t = [1:5*60:2*3600*24]; %Time vector in s - take position every 5 mins
t = linspace(0,10*day_sd, 100);%linearly spaced tvec points in s (better for analysis, no periodic behavior). 


%NOTE: For retrograde burns, thrust should be negative
max_thrust = 330e-6; %nominal thrust for Enpulsion Nano (https://www.enpulsion.com/wp-content/uploads/ENP2018-001.G-ENPULSION-NANO-Product-Overview.pdf)
I_sp = 2000; %specific impulse in s
prop_power = 40;% Power of electric propulsion system in W
m_0 = 8; %spacecraft wet mass in kg
m_dot = max_thrust/(9.81*I_sp); % mass flow rate from prop %%!!! Modify for actual thrust in propellant consumtion

r_0 = a; %initial orbital radius
r_f = 497.5 + Re;
alt_tol = 0.2; %Tolerance for final altitude [km]. Default is 0.2 km. 
% t_ltt = (m_0/(exp(((1/sqrt(r_0))-(1/sqrt(r_f)))*(sqrt(mu)*m_dot/thrust)))-m_0)/(-m_dot); %Time of the low thrust transfer to achieve orbit transition from desired GOM to ROM
% !!! CHECK THE ABOVE LINE!!! GET NEGTIVE NUMBER HERE (WRONG)


[posN,RAAN, RAAN_dot ,v,alt,thrust, t] = PropagateOrbit_J2(t,mu,a,e,h,incl,RAAN,ArgPer,anomaly, max_thrust, m_0,m_dot,prop_power, r_f, alt_tol);


figure()
%make and plot sphere for earth
[sx,sy,sz] = sphere();
r = 6378;
surf( r*sx, r*sy, r*sz )
hold on
axis equal
%plot propagated orbit at times t
plot3(posN(1,:),posN(2,:),posN(3,:),'r.');

xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
hold off

figure()
plot(a, rad2deg(omega_E)./rad2deg(RAAN_dot));
xlabel('Semi-major Axis [km]')
ylabel('Omega_E/RAAN dot');
%ideally, we want the rotation rate of the earth to be an integer multiple of
%the rate of precession of the orbit
hold on
hold off

% for i = 2:length(RAAN_dot)
%     trapz(t(1:i), RAAN_dot(1:i));
% end

figure()
plot(t, RAAN)
xlabel('Time [s]')
ylabel('Orbital Precession [deg]')

figure()
plot(t, a)
xlabel('Time (s)')
ylabel('Semi-major Axis [km]');

% Need special toolbox to run the following
figure()
fnplt(cscvn(posN)); hold on
plot3(posN(1,:), posN(2,:), posN(3,:), 'o')



