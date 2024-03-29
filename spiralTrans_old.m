%% Function Description
% 
% Perform analytical approximations for spiral transfer through numerical
% iteration. Use PropagateOrbit_J2 for more in-depth solutions with orbital
% position calculated for each point in time. 
%
% William Parker - March 2021

function [t_trans, dv, num_trans] = spiralTrans_old(initial_alt, dv_avail,e,incl, delta_RAAN, Re, mu, J2, dv_per_trans, accel)

num_days = 40;
iter = num_days*24; %number of iterations (1 per hour)
dt = 3600; %dt is one hour
tstep = [0,dt,iter*dt];


r0 = initial_alt + Re;
v0 = sqrt(mu/r0);

RAAN_dot_init = -((3/2)*(sqrt(mu)*J2*Re^2)/((1-e^2)^2*r0^(7/2)))*cosd(incl); %Nodal precession in deg/s for initial orbit


if delta_RAAN > 0
    a = accel; 
else 
    a = -accel;
end

r = zeros(iter, 1);
RAAN_dot = zeros(iter,1);
RAAN_drift = zeros(iter, 1);
t = zeros(iter,1);

for i = 1:iter
    
    if RAAN_drift(i) <= abs(0.5*delta_RAAN)
        a = a; %acceleration is in positive direction
    elseif RAAN_drift(i) > abs(0.5*delta_RAAN) && RAAN_drift(i) <= delta_RAAN
        a = -1*a; %acceleration is in the negative direction
    else
        a = 0; %don't apply any thrust
        t_trans = t;
        dv = 2*abs(sqrt(mu/(max(r)))-sqrt(mu/min(r)));
        num_trans = dv_avail/dv;
        return
    end
    
    r(i) = r0/(1-a*dt/v0)^2;
    RAAN_dot(i) = -((3/2)*(sqrt(mu)*J2*Re^2)/((1-e^2)^2*r(i)^(7/2)))*cosd(incl); %Nodal precession in deg/s for initial orbit

    if i > 1
        RAAN_drift(i) = RAAN_drift(i-1)+ (RAAN_dot(i)-RAAN_dot_init)*dt; 
        t(i) = t(i-1)+ dt; 
    end


    r0 = r(i);
    v0 = sqrt(mu/r0);
    
end
figure()
plot(t, r)

figure()
plot(t, RAAN_drift)

end

    
    
    