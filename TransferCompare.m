%% Function Description
% 
% Run this function from Run_TransferCompare.m. In that script, you can
% iterate between transfer types and the desired change in RAAN to help
% determine which transfer type is best for a given new target location. 
%
% Developed by William Parker - March 2021

function [t_trans, dv, num_trans] = TransferCompare(trans_type, thrust_level, delta_RAAN_or, idx)

[planet,t,orbInit,orbFin,numHT,numLT,sc,prop] = AssignParams(thrust_level);

if idx > 2 %only do this the first iteration 
    save('planet.mat', '-struct', 'planet');
    save('orbInit.mat', '-struct', 'orbInit');
    save('orbFin.mat', '-struct', 'orbFin');
    save('sc.mat', '-struct', 'sc');
    save('prop.mat', '-struct', 'prop');
end

    load('planet.mat');
    load('orbInit.mat');
    load('orbFin.mat');
    load('sc.mat');
    load('prop.mat');

delta_RAAN = delta_RAAN_or;


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
    posinfo = 2; %don't want to propagate orbit position (takes a long time, and isn't useful for comparison).
    [posN,RAAN, RAAN_dot, RAAN_rel, ArgPer,v,alt,thrust, t, t_trans, dv, num_trans] = spiralTrans(t,mu,a,e,h,incl,RAAN,ArgPer,anomaly, thrust, m_0,m_dot,prop_power, r_f, alt_tol, RAAN_f, Re, dv_avail, posinfo);
    fprintf('<strong>Time-optimal spiral transfer completed (using low thrust)!</strong>\n')

elseif trans_type == 4 %Direct Plane Change
    [t_trans, dv, num_trans] = PlaneChangeTrans(initial_alt, dv_avail,e,incl, delta_RAAN, Re, mu, J2);
    fprintf('<strong>Direct plane change transfer completed (using high thrust)!</strong>\n')

elseif trans_type == 5 % Low thrust Plane Change
    [t_trans, dv, num_trans] = planeChange_LT(initial_alt, dv_avail,incl, delta_RAAN, Re, mu, accel);
    fprintf('<strong>Low thrust plane change transfer completed!</strong>\n')

end

fprintf('<strong>Transfer time:</strong> %.2f seconds (or %.3f hours or %.3f days) \n', t_trans, t_trans/3600, t_trans/(day_sd));
fprintf('<strong>Delta-v required:</strong> %.4f km/s or %.2f m/s \n', dv, dv*1000);
fprintf('<strong>Number of transfers possible with this prop system and transfer type:</strong> %.3f \n', num_trans);



end


