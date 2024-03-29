%% WELCOME TO THE ROAMS SIMULATION TOOL (RST AKA "RUSTY")
% 
% Use this script to run ROAMS maneuver simulations. First, input planetary properties, 
% simulation timing, current & desired orbit elements, spacecraft
% properties, and propulsion system parameters in AssignParams.m. The, run
% this script to generate a simulation for various methods of orbital transfer, 
% and produce performance tracking metrics in the process for analysis. To
% compare multiple transfer methods over a variety of transfer distances,
% use Run_TransferCompare.m.
% 
% DIRECTIONS FOR USE: 
% 1. Edit the input parameters in AssignParams.m for desired spacecraft/orbit
% 2. Answer any questions asked in input dialog boxes or in the command
% window (HT/LT = High/Low Thrust) -- These may not pop up automatically,
% so if the simulation is paused for an extended period of time, look in
% your matlab tabs for a selection menu window. 
% 3. Wait for orbit propagator to work through the entire simulation (it
% will continue for the simulation duration you assigned, and not stop
% automatically when the transfer has been completed). 
% 4. View output information for selected transfer type (should be fast to compute for
% most transfer types)
% 5. NOTE: If you're simulating a spiral transfer and get an error that
% says "Increase the duration of your simulation window," increase numDays in AssignParams.m.
% If your spiral transfer plots don't appear to perform as expected, try 
% increasing numPoints. The spiral transfer simulation is especially
% susceptible to quantization error.
% 6. If you would like to compare different transfer methods over multiple
% transfer scenarios (distance between targets is analagous to delta_RAAN),
% edit the parameters in AssignParams.m and then run
% Run_TransferCompare.m.
% 
% 
% ASSUMPTIONS: 
% 1. Only circular orbits in spiral transfers
% 2. No inclination change in transfer, only altitude changes (except for direct plane change maneuver). 
% 3. Constant mass 
% 4. Perfect burns - exactly prograde/retrograde at nominal thrust
% 5. No environmental disturbance forces (drag, SRP, etc.)
% 6. Impulsive (instantaneous) burns for high thrust transfers.
% 7. No RAAN drift during high-thrust transfers to/from GOM (only while in GOM). 
% 
% CURRENT LIMITATIONS:
% 1. Orbit propagator is only providing data for low thrust transfers. High
% thrust transfers are simulated using analytical approximations, rather
% than numerical iteration through time (for impulsive maneuvers, the
% maneuver itself is assumed to occur over a very small timescale)
% 2. Assuming purely impulsive high thrust transfers
%
% Developed by William Parker - March 2021

%% 
clc; clear all; close all;
%% PROMPT USER FOR SIMULATION TYPE
% HT = High Thrust, LT = Low Thrust
[trans_type] = menu('What type of orbital transfer would you like to simulate?','Circular GOM Transfer (HT)', 'Elliptical Transfer (HT)', 'Optimized Spiral Transfer (LT)', 'Plane Change (HT)', 'Plane Change (LT)');

if trans_type == 1 
    thrust_level = 1;
    disp('For a circular GOM transfer, assume high thrust impulsive maneuvers. For a low thrust system, try the spiral transfer.');
end

if trans_type == 2 || trans_type == 4
    thrust_level = 1;
end

if trans_type == 3 || trans_type == 5
    thrust_level = 2;
end

%% ASSIGN PARAMETERS
[planet,t,orbInit,orbFin,numHT,numLT,sc,prop] = AssignParams(thrust_level);

save('planet.mat', '-struct', 'planet');
save('orbInit.mat', '-struct', 'orbInit');
save('orbFin.mat', '-struct', 'orbFin');
save('sc.mat', '-struct', 'sc');
save('prop.mat', '-struct', 'prop');

load('planet.mat');
load('orbInit.mat');
load('orbFin.mat');
load('sc.mat');
load('prop.mat');

%% PROPAGATE SCENARIO WITH J2 CORRECTION
%For low thrust transfers, there is also an option to propagate the orbit
%and determine the spacecraft position over time. This is not necessary for
%most design tasks and takes significantly longer to perform the
%simulation. 

if trans_type == 1 %GOM transfer
    [t_trans, dv, num_trans] = GOMTrans(GOM_alt_diff, initial_alt, final_alt, dv_avail, e, incl, delta_RAAN, Re, mu, J2);
    fprintf('<strong>Circular GOM transfer completed (using high thrust)!</strong>\n')

elseif trans_type == 2 %Elliptical transfer
    %assign dv we're willing to sacrifice per maneuver (>dv, <t)
    dv_per_trans = dv_avail/numHT; %divide total dv by number of transfers we'd like
    [t_trans, dv, num_trans] = ellipticalTrans(initial_alt, dv_avail,e,incl, delta_RAAN, Re, mu, J2, dv_per_trans);
    fprintf('<strong>Elliptical transfer completed (using high thrust)!</strong>\n')

elseif trans_type == 3 %low thrust spiral transfer
    dv_per_trans = dv_avail/numLT;
    posinfo = menu('Propagate spacecraft position (significantly longer simulation time)?', 'Yes', 'No');
    [posN,RAAN, RAAN_dot, RAAN_rel, ArgPer,v,alt,thrust, t, t_trans, dv, num_trans] = spiralTrans(t,mu,a,e,h,incl,RAAN,ArgPer,anomaly, thrust, m_0,m_dot,prop_power, r_f, alt_tol, RAAN_f, Re, dv_avail, posinfo);
    fprintf('<strong>Time-optimal spiral transfer completed (using low thrust)!</strong>\n')

elseif trans_type == 4 %High thrust Plane Change
    [t_trans, dv, num_trans] = PlaneChangeTrans(initial_alt, dv_avail,e,incl, delta_RAAN, Re, mu, J2);
    fprintf('<strong>High thrust plane change transfer completed!</strong>\n')

elseif trans_type == 5 % Low thrust Plane Change
    [t_trans, dv, num_trans] = planeChange_LT(initial_alt, dv_avail,incl, delta_RAAN, Re, mu, accel);
    fprintf('<strong>Low thrust plane change transfer completed!</strong>\n')
end
fprintf('<strong>Change in RAAN: </strong> %.3f deg \n', delta_RAAN);
fprintf('<strong>Transfer time:</strong> %.2f seconds (or %.3f hours or %.3f days) \n', t_trans, t_trans/3600, t_trans/(day_sd));
fprintf('<strong>Delta-v required:</strong> %.4f km/s or %.2f m/s \n', dv, dv*1000);
fprintf('<strong>Number of transfers possible with this prop system and transfer type:</strong> %.3f \n', num_trans);


%if low thrust spiral transfer, generate results from simulation iterations
if trans_type == 3 
    %% RECORD ORBITAL ELEMENTS AT EACH TIMESTEP OF SIMULATION
    % Format [e, a, i, RAAN, ArgPer]
    % Note: anomaly varies over an orbit, just need these 5 elements to define the
    % orbit shape and orientation
    orbElem(:,1) = zeros(length(t),1);
    orbElem(:,2) = alt'+Re;
    orbElem(:,3) = ones(length(t),1).*incl;
    orbElem(:,4) = rad2deg(RAAN');
    orbElem(:,5) = rad2deg(ArgPer);

    %% PLOT RESULTS
    % convert time from s to days
    t = t/day_sd;
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
%         figure()
%         plot(alt, rad2deg(omega_E)./rad2deg(RAAN_dot));
%         xlabel('Altitude [km]')
%         ylabel('Omega_E/RAAN dot');
%         %ideally, we want the rotation rate of the earth to be an integer multiple of
%         %the rate of precession of the orbit
%         hold off

    %Plot RAAN of orbit over time
        figure()
        plot(t, rad2deg(RAAN))
        xlabel('Time [day]')
        ylabel('Orbit Longitudinal Precession(RAAN)[deg]')

    %Plot altitude over time
        figure()
        plot(t, alt)
        xlabel('Time [day]')
        ylabel('Altitude [km]');

    %Plot Thrust over time
        figure()
        plot(t, thrust)
        xlabel('Time [day]')
        ylabel('Commanded Thrust [N]')

    % Plot relative RAAN over time (RAAN precession from phasing maneuver)
        figure()
        plot(t, rad2deg(RAAN_rel))
        xlabel('Time [day]');
        ylabel('Relative RAAN [deg]')
        yline(delta_RAAN, 'r');
        legend('Relative RAAN change', 'Commanded RAAN change');
end



