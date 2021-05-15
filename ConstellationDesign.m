%% Compare different constellation architectures and the performance requirements for the spacecraft
% William Parker, May 2021
% clc;clear all; close all; 
%% Set Requirements
load('planet.mat');
% Req_t2targ = [1:0.1:25].*86164.1; %Time between targets in seconds
Req_t2targ = [7].*86164.1; %Time between targets in seconds
Req_numTrans = 20; %Number of transfers that are possible for each spacecraft (worst case scenario)
Req_tRevisit = 3*3600; %revisit time in seconds

%% Analysis for number of orbital planes
% For sunsynch (561 km, 97.64 km), ~86200s before repeat ground track
% (almost exactly one day)
% @ equator (max case), ~24 deg separation between ground tracks (max
% distance to travel is +/- 12 deg if spaced evenly). 

% time to target requirement translates directly to number of orbital
% planes (since spacing between planes makes transitions take
% longer/shorter). 

numPlanes = [8];
dist_per_sat = zeros(length(numPlanes),1); %max angular RAAN distance that a satellite needs to travel in deg 
RAAN_drift_rate = zeros(length(numPlanes),1);  % drift rate in deg/s
for i = 1:length(Req_t2targ)
    dist_per_sat(i) = 12./numPlanes; %in deg
    RAAN_drift_rate(i) = dist_per_sat(i)./Req_t2targ(i);
end
% subplot(1,3,1)
% plot(numPlanes,RAAN_drift_rate*86164.1, 'k-o')
% xlabel('Number of Planes');
% ylabel('Required RAAN Drift Rate [deg/day]');
% hold on

% For high thrust, assume hohmann transfer. Calculate ROM-GOM altitude
% change required for circular GOM to achieve the required drift rate. 
a_init = 497+Re; % orbit radius in km
% incl = 97.64; %inclination in degrees
incl = 51.6;
e = 0;

RAAN_dot_init = -((3/2)*(sqrt(mu)*J2*Re^2)/((1-e^2)^2*a_init^(7/2)))*cosd(incl); %Nodal precession in rad/s for initial orbit
RAAN_dot_fin = RAAN_dot_init+deg2rad(RAAN_drift_rate);
% From final RAAN_dot, solve for final altitude (to get ROM-GOM alt
% difference)
% a_fin = ((-2.*RAAN_dot_fin.*1)./(3.*sqrt(mu).*J2.*Re.^2.*cosd(incl))).^(-2/7);
a_fin = ((-2/3).*RAAN_dot_fin./(cosd(incl)*sqrt(mu)*J2*Re.^2)).^(-2/7);
alt_change = abs(a_fin-a_init);

% subplot(1,3,2)
% plot(numPlanes, alt_change, 'x-');
% xlabel('Number of Planes');
% ylabel('Necessary High Thrust Altitude Change [km]');
% hold on

% Calculate delta-v requirements for high thrust onboard in order to
% achieve 20 transfers to the altitude required for the prescribed number
% of orbital planes

% Assume hohmann transfer
dv1 = sqrt(mu./a_init).*(sqrt(2.*a_fin./(a_init+a_fin))-1);
dv2 = sqrt(mu./a_fin).*(1-sqrt(2*a_init./(a_init+a_fin)));
dv_outNback = abs(2.*(dv1+dv2))*1000;

% subplot(1,3,3)
% plot(numPlanes, dv_outNback,'*-');
% xlabel('Number of Planes');
% ylabel('Delta-V Required Per Transfer [m/s]');
% hold on

%% For Elliptical Transfer
%can change both a and e
alt_diff = [1:1:200];
rp = a_init;
ra = rp+alt_diff;
a_fin = (rp+ra)./2;
e = (ra-rp)./(ra+rp);
RAAN_dot_ellip = -((3/2).*(sqrt(mu)*J2*Re^2)./((1-e.^2).^2.*a_fin.^(7/2))).*cosd(incl);

% plot(alt_diff, RAAN_dot_ellip)
% hold on
% yline(RAAN_dot_fin, 'g');
% yline(RAAN_dot_init, 'r');



%% FOR LOW THRUST ANALYSIS
% Want to determine delta-v, thrust level, and duration required per transfer

% Use triangle geometry to define max difference in RAAN for low thrust
% (approximate to be linear over altitude changes) A_tri = bh/2
max_relRAAN = 2.*dist_per_sat./Req_t2targ; % in degrees
max_relRAAN_rad = deg2rad(max_relRAAN); % in rad

% determine altitude of max_relRAAN
RAAN_dot_init = -((3/2)*(sqrt(mu)*J2*Re^2)/((1-e^2)^2*a_init^(7/2)))*cosd(incl); %Nodal precession in rad/s for initial orbit
RAAN_dot_fin = RAAN_dot_init-max_relRAAN_rad;
a_fin = ((-2/3).*RAAN_dot_fin./(cosd(incl)*sqrt(mu)*J2*Re.^2)).^(-2/7);
alt_change = abs(a_fin-a_init);

% Determine the thrust it takes to achieve this altitude change in the
% prescribed time requirement with a low thrust transfer
mass = 12; %kg
accel = -(sqrt(mu./a_init)./Req_t2targ).*((a_init./(a_init+alt_change)).^2 - 1);
thrust_LT = accel*(1000)*mass; %1000 converts km to m, so thrust is in N
dV_LT = accel.*1000.*Req_t2targ;

subplot(2,1,1)
semilogy(Req_t2targ./86164.1, thrust_LT, '-');
xlabel('Requirment for Reconfigure Time [days]');
ylabel('LT Thrust Required for Transfer Time Req. [N]');
hold on

subplot(2,1,2)
semilogy(Req_t2targ./86164.1, dV_LT, '-');
xlabel('Requirement for Reconfigure Time [days]');
ylabel('LT \DeltaV Per Maneuver');
hold on







