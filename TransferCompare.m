function [t_trans, dv, num_trans] = TransferCompare(trans_type, thrust_level, delta_RAAN)

% %% PROMPT USER FOR SIMULATION TYPE
% % HT = High Thrust, LT = Low Thrust
% [trans_type] = menu('What type of orbital transfer would you like to simulate?','Circular GOM Transfer (HT/LT)', 'Elliptical Transfer(HT)', 'Optimized Spiral Transfer (LT)', 'Direct Plane Change (HT)');
% 
% if trans_type == 1
%     [thrust_level] = menu('High or Low Thrust? Be sure to have entered thruster properties and GOM altitude in the main script.', 'High Thrust', 'Low Thrust');
% end
% 
% if trans_type == 2 || trans_type == 4
%     thrust_level = 1;
% end
% 
% if trans_type == 3
%     thrust_level = 2;
% end

%% SPECIFY PLANETARY PROPERTIES
Re = 6378.15; %Radius of Earth in km
omega_E = 7.292e-5; %rotation rate of the earth in rad/s
day_sd = 86164.1; %number of seconds in a sidereal day
mu = 398600; %Grav. param. for Earth in km^3/s^2
g0 = 9.81; % Gravitational acceleration in m/s^2;
J2 = 1.08263e-3; 

%% SPECIFY SIMULATION TIME VECTOR
% t = [1:5*60:2*3600*24]; %Time vector in s - take position every 5 mins
t = linspace(0,10*day_sd, 100);%linearly spaced tvec points in s (better for analysis w/o plotting, no periodic behavior). 

%% SPECIFY INITIAL ORBIT ELEMENTS
%ROAMS GOM orbit parameters (from 16.851 Fall 2021 Final Pres.)
initial_alt = 497.5; %initial altitude in km
e = 0;
incl = 51.6; % inclination in deg
RAAN = 60; % in deg. initially say 0 for placeholder
ArgPer = 0; % in deg. initially say 0 for placeholder
anomaly = 0; % in deg. initially say 0 for placeholder

a = initial_alt+Re;
h = sqrt(a*mu);
T = 2*pi/(sqrt(mu))*a^(3/2); %Orbital period in s

%% SPECIFY FINAL DESIRED ORBIT ELEMENTS
final_alt = 497.5; %in km, should be same as initial_alt for ROM->ROM transfers
RAAN_f = 70; %in deg. 
% delta_RAAN = RAAN_f-RAAN; %difference in deg. NOTE: Revisit for changeover back to 0 (prevent differencing issue). 
alt_tol = 0.2; %Tolerance for final altitude [km]. Default is 0.2 km. %NOTE: May not see this work well for large timesteps

GOM_alt_diff = 50; % ROM-GOM altitude difference in km (if simulating GOM transfers). We can go +diff or -diff to move in different RAAN directions.

r_0 = a; %initial orbital radius
r_f = final_alt + Re;

%% SPECIFY CONSTRAINTS ON MANEUVERS
% t_f = 10;%Time in DAYS
% t_f = t_f*day_sd; %convert to seconds

numHT = 5; %number of desired ROM->ROM maneuvers with high thrust
numLT = 2; %number of desired ROM->ROM maneuvers with low thrust

%% SPECIFY SPACECRAFT PARAMETERS
m_0 = 12; %spacecraft wet mass in kg

%% SPECIFY PROPULSION SYSTEM PARAMETERS
if thrust_level == 1
    %HIGH THRUST
    thrust = 1; %nominal thrust for Aeroject GR1 in N 
    I_sp = 220; % (May not be the correct number for Aerojet GR1)
    prop_mass = 0.5; % in kg for 1U Aerojet GR1
    prop_power = 20; %Power in W (assumed 20 - not necessarily accurate for GR1)
elseif thrust_level == 2
    %LOW THRUST
    thrust = 330e-6; %nominal thrust for Enpulsion Nano in N (https://www.enpulsion.com/wp-content/uploads/ENP2018-001.G-ENPULSION-NANO-Product-Overview.pdf)
    I_sp = 6000; %specific impulse in s for impulsion Nano (high estimate)
    prop_power = 40;% Power of electric propulsion system in W
    prop_mass = 0.220; % initial mass of propellant in kg
end

accel = thrust/m_0;
m_dot = thrust/(I_sp*g0); % mass flow rate in kg/s;
t_thrust = prop_mass/m_dot; % total time in s that thruster can burn before running out of propellant
dv_avail = t_thrust*accel/1000; %dv available in km/s

% Note: For the current implementation, we assume impulsive maneuvers and
% constant mass, so specifics for the high thrust system are only used to
% determine maximum delta-v. 

%% PROPAGATE SCENARIO WITH J2 CORRECTION
% At the moment, we're only able to propagate low thrust transfers (high
% thrust are done from analytical approximations, not numerical
% simulation).

%if low thrust... propagate orbit and generate plots
if thrust_level == 2
   [posN,RAAN, RAAN_dot, RAAN_rel, ArgPer,v,alt,thrust, t] = PropagateOrbit_J2(t,mu,a,e,h,incl,RAAN,ArgPer,anomaly, thrust, m_0,m_dot,prop_power, r_f, alt_tol, RAAN_f, Re);
end

%if high thrust... analytically determine transfer time, dv requirement,
%number of transfers possible
if thrust_level == 1
    if trans_type == 1 %GOM transfer
        [t_trans, dv, num_trans] = GOMTrans(GOM_alt_diff, initial_alt, final_alt, dv_avail, e, incl, delta_RAAN, Re, mu, J2);
        fprintf('<strong>Circular GOM transfer completed!</strong>\n')
        
    elseif trans_type == 2 %Elliptical transfer
        %assign dv we're willing to sacrifice per maneuver (>dv, <t)
        dv_per_trans = dv_avail/numHT; %divide total dv by number of transfers we'd like
        [t_trans, dv, num_trans] = ellipticalTrans(initial_alt, dv_avail,e,incl, delta_RAAN, Re, mu, J2, dv_per_trans);
        fprintf('<strong>Elliptical transfer completed!</strong>\n')
        
    elseif trans_type == 4 %Direct Plane Change
        [t_trans, dv, num_trans] = PlaneChangeTrans(initial_alt, dv_avail,e,incl, delta_RAAN, Re, mu, J2);
        fprintf('<strong>Direct plane change transfer completed!</strong>\n')
    end
    
    fprintf('<strong>Transfer time:</strong> %.2f seconds (or %.3f hours or %.3f days) \n', t_trans, t_trans/3600, t_trans/(day_sd));
    fprintf('<strong>Delta-v required:</strong> %.4f km/s or %.2f m/s \n', dv, dv*1000);
    fprintf('<strong>Number of transfers possible with this prop system and transfer type:</strong> %.3f \n', num_trans);
end


if thrust_level == 2
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
        plot(t, rad2deg(RAAN))
        xlabel('Time [s]')
        ylabel('Orbit Longitudinal Precession(RAAN)[deg]')

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

    % Plot relative RAAN over time (RAAN precession from phasing maneuver)
        figure()
        plot(t, rad2deg(RAAN_rel))
        xlabel('Time [s]');
        ylabel('Relative RAAN [deg]')
end
end


