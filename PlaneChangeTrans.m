%% Function Description
% Determine the performance metrics for changing orbital RAAN directly in
% ROM using a plane change maneuver
% 
% William Parker - March 2021

function [t_trans, dv, num_trans] = PlaneChangeTrans(initial_alt, dv_avail,e,incl, delta_RAAN, Re, mu, J2)

v = sqrt(mu/(initial_alt+Re));

%use law of sins to determine velocity change requirement
ang = (180-(delta_RAAN))/2;
dv = sind(delta_RAAN)*v/sind(ang);

t_trans = 0; %we're going assume impulsive maneuvers here (zero time to carry out maneuver)

num_trans = dv_avail/dv;
end
