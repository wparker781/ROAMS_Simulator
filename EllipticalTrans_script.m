clc; clear all; close all; 
%% Script for identifying optimal trajectory to achieve a given orbital transition
% Objective function: time, propellant, delta-v, power consumed
% Constraints: desired change in RAAN (given)
mu = 398600;
J2 = 1.08263e-3; 
Re = 6378.15;
% e = [0,0.01,0.1];
e = 0;
% alt = [497.5:10:550];
alt = 497.3;
a = alt+Re;
incl = 51.6;%deg

omega_E = 7.292e-5; %rotation rate of the earth in rad/s

RAAN_dot_ROM = -((3/2).*(sqrt(mu).*J2.*Re.^2)./((1-e.^2).^2.*a.^(7/2))).*cos(deg2rad(incl)); %Nodal precession in rad/s


% n = [1:1:100];
% RAAN_dot = -n*omega_E + omega_E; 
% % t1 = -2.*RAAN_dot.*cos(deg2rad(incl)).*1;
% % t2 = (3.*sqrt(mu).*J2.*Re.^2);
% % a = (t1./t2)
% % a = a.^(-(2/7))
% a = ((-2.*RAAN_dot.*cos(deg2rad(incl)).*1)./(3.*sqrt(mu).*J2.*Re.^2)).^(-2/7);
% alt = a-Re; 
% plot(n, alt)

%For an eliptical impulsive transfer, just do first half of a hohmann transfer
%(and double it for the delta-v required to return to the circular orbit). 

% Assign the delta-v that we're willing to sacrifice for one transition
% (ROM->ROM)
dv_per_trans = 0.090; % in KM/s
% Increasing altitude or decreasing altitude? increase_alt = 1 for yes, 0
% for no. 
increase_alt = 0; 
dv_diff = dv_per_trans/2; %diference in v between ROM and transfer orbit
if increase_alt == 1 
    %maneuver will occur at transfer orbit periapsis
    rp = a;
    v_per = dv_diff + sqrt(mu/a);
    ra = (v_per^2*rp^2)/(2*mu-v_per^2*rp);
    e = (ra-rp)/(ra+rp);
    a = (ra+rp)/2;
    RAAN_dot_trans = -((3/2).*(sqrt(mu).*J2.*Re.^2)./((1-e.^2).^2.*a.^(7/2))).*cos(deg2rad(incl)); %Nodal precession in rad/s
            
    else if increase_alt == 0
        %maneuver will occur at transfer orbit apoapsis
        ra = a;
        v_apo = sqrt(mu/a)-dv_diff;
        rp = v_apo^2*ra^2/(2*mu-v_apo^2*ra);
        e = (ra-rp)/(ra+rp);
        a = (ra+rp)/2;
        RAAN_dot_trans = -((3/2).*(sqrt(mu).*J2.*Re.^2)./((1-e.^2).^2.*a.^(7/2))).*cos(deg2rad(incl)); %Nodal precession in rad/s
        
        else 
            disp('increase_alt may only be assigned 0 or 1'); 
        end
end

            
%Now, find difference in RAAN_dot between transfer orbit and ROM. Generate
%drift rate in deg/day;

RAAN_drift_rad = RAAN_dot_trans-RAAN_dot_ROM; %in Rad/s
RAAN_drift_deg = rad2deg(RAAN_drift_rad);%in deg/s
degPerDay = RAAN_drift_deg*86164.1 % deg/day




RAAN_dot_rad = -((3/2).*(sqrt(mu).*J2.*Re.^2)./((1-e.^2).^2.*a.^(7/2))).*cos(deg2rad(incl)); %Nodal precession in rad/s
nodal_T = 2*pi./(RAAN_dot_rad);
check = nodal_T/((2*pi)/omega_E); 
n = (RAAN_dot_rad-omega_E)/omega_E
RAAN_dot_deg = rad2deg(RAAN_dot_rad); %in degrees/s


RAAN_dot_hour = RAAN_dot_deg*3600; %Nodal precession in deg per hour

plot(e, RAAN_dot_hour)

function [vp] = vp(mu,ra,rp)
    vp = sqrt((2*mu*ra)/(rp*(rp+ra)));
end

function [va] = va(mu, ra,rp)
    va = sqrt((2*mu*rp)/(ra*(rp+ra)));
end

function ra = r_apo(mu, v_per, rp)
    syms ra
    solve(v_per == sqrt((2*mu*rp)/(ra*(rp+ra))));
    ra = double(ra);
end



