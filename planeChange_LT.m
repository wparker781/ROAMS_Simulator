%% Function Description
% Perform direct plane change maneuver using low thrust. See 
% https://apps.dtic.mil/sti/pdfs/ADA384536.pdf, pg. 19 for equation and
% derivation. 
%
% Developed by William Parker - March 2021
%
function [t_trans, dv, num_trans] = planeChange_LT(initial_alt, dv_avail,incl, delta_RAAN, Re, mu, accel)

dv = (pi/2)*(sqrt(mu/(initial_alt+Re)))*sin(deg2rad(incl))*abs(deg2rad(delta_RAAN));
num_trans = dv_avail/dv;
%to determine the time it takes to carry out this maneuver, just consider
%that we know the dv and we know the acceleration rate of the spacecraft.
%Determine the time it takes for the spacecraft to accumulate that change
%in velocity at the given acceleration rate (assume constant thrusting
%throughout maneuver). 
t_trans = dv/accel;