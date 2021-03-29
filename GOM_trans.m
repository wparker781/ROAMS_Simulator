%% FUNCTION DESCRIPTION
% Perform a Hohmann transfer from the initial RGT orbit to GOM and back.
% Determine the total delta-V required to perform this manneuver, and the
% number of transfers that would be possible with the proposed high thrust
% propulsion system. Also determine the time it would take to perform the
% transfer (including the homann transfer and the time it takes for the
% RAAN precession due to J2 to accumulate to the desired angular distance).
% For our purposes, assume that there is no RAAN drift during the hohmann
% transfer (should be small because transfer is fast anyways). 
% 
% William Parker - March 2021

function [t_tot, dv, num_trans] = GOMTrans(GOM_alt_diff, initial_alt, final_alt, dv_avail, e, incl, delta_RAAN, Re)
%check to see which direction we want to move in for RAAN change (GOM
%should be +/- depending on which direction we want to move in)
if delta_RAAN >= 0
    GOM_alt = initial_alt + GOM_alt_diff; 
else
    GOM_alt = initial_alt - GOM_alt_diff;
end

% perform hohmann transfer - figure out delta v and duration for transfer
r1 = initial_alt + Re;
r2 = GOM_alt + Re;
Vpto = sqrt(mu*((2/(r1))-(1/(r2))));
Vato = Vpto*(r1/r2);
dv1 = Vpto-sqrt(mu/r1);
dv2 = sqrt(mu/r2)-Vato;
dv_out = dv1 + dv2; %this is delta v (in m) one way. For 2-way, we need to add the way back.

r1 = GOM_alt + Re;
r2 = final_alt + Re;
Vpto = sqrt(mu*((2/(r1))-(1/(r2))));
Vato = Vpto*(r1/r2);
dv1 = Vpto-sqrt(mu/r1);
dv2 = sqrt(mu/r2)-Vato;
dv_back = dv1 + dv2;

dv = dv_out + dv_back; %dv in km/s

num_trans = dv_avail/dv;

%For hohmann transfer time, just consider elliptical orbit period (half
%out, half back). 
a = r1+r2;
t_ht = 2*pi*a^(3/2)./sqrt(mu); %hohmann transfer time (both ways) in s

%Now, determine the amount of time required for the RAAN drift to reach
%the desired magnitude
a = initial_alt;
RAAN_dot_init = -((3/2)*(sqrt(mu)*J2*Re^2)/((1-e^2)^2*a^(7/2)))*cosd(incl); %Nodal precession in deg/s for initial orbit
a = GOM_alt;
RAAN_dot_GOM = -((3/2)*(sqrt(mu)*J2*Re^2)/((1-e^2)^2*a^(7/2)))*cosd(incl); %Nodal precession in deg/s for GOM orbit

RAAN_drift_rate = RAAN_dot_GOM-RAAN_dot_init;
t_drift = delta_RAAN/RAAN_drift_rate; %drift time in s

t_tot = t_ht+t_drift;

end

