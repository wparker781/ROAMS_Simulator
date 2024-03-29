function [posN,RAAN_track,RAAN_dot_track,RAAN_rel_track,ArgPer_track, v_track, alt_track, thrust_track, t] = spiralTrans(t,mu,a,e,h,incl,RAAN,ArgPer,anomaly, thrust_max, sc_mass, m_dot,prop_power, r_f, alt_tol, RAAN_f, Re)
%% Function Description
% Perform low thrust spiral transfer for time-optimal maneuver. This
% function runs iteratively to track relative RAAN, altitude, etc over
% time. 
% 
% INPUTS: 
% t = nx1 array of times since initial actual anomaly (s)
% mu = gravitational parameter of the attracting body
% a = semimajor axis 
% e = eccentricity of orbit
% h = angular momentum
% incl = inclination (degrees)
% RAAN = right ascention of the ascending node (degrees)
% ArgPer = argument of periapsis (degrees)
% anomaly = actual anomally at the first ephemeris time of measurement(deg)
% thrust = thrust from electric propulsion system for low-thrust
%   maneuvering. For this orbit propagation tool, we assume continuous thrust
%   over the duration of the input time vector (N). 
% sc_mass = wet mass of spacecraft
% r_f = semi-major axis of final desired orbit

% OUTPUTS:
% PosN = an nx3 matrix of position vectors with 3 components for each time t [km]
% RAAN (updated at each timestep because of the J2 perturbation)
% t = nx1 array of times since initial actual anomaly (s)

% NOTE: Assuming Earth for radius and J2
% 
% Developed by William Parker, November 2019
% 3/1/2021 - Modified to add J2 perturbation and ROAMS-relevant altitude
% controller
% 3/29/2021 - Modified to gear towards low thrust spiral transfers
% specifically

%% Initialize variables

% solE = zeros(length(t),1);
% theta = zeros(length(t),1);
% r = zeros(length(t),1);
% PosOrb = zeros(3,length(t));
% posPF = zeros(3,length(t));
posN = zeros(3,length(t));
RAAN_track = zeros(1,length(t));
RAAN_dot_track = zeros(1,length(t)); 
ArgPer_track = zeros(1,length(t));
v_track = zeros(1,length(t)); 
alt_track = zeros(1,length(t));
thrust_track = zeros(1,length(t)); 
delta_v_track = zeros(1,length(t)); 
flag = 0; 
resDisp = 0; 

J2 = 1.08263e-3; 
dt = diff(t);


%convert angles in degrees to radians
incl = deg2rad(incl);
RAAN = deg2rad(RAAN);
ArgPer = deg2rad(ArgPer);
anomaly = deg2rad(anomaly);

%Calculate velocity  of sc in the initial orbit
v = sqrt(mu/a);
alt_f = r_f-Re;
alt_0 = a-Re;
int_err = 0;
prev_err = r_f - a;
prev_alt = alt_0;
thrust = 0; %for the initial iteration, don't use any thrust. Just measure error. 
delta_v = 0; %no delta-v on the first iteration
RAAN_0 = RAAN;
RAAN_f = deg2rad(RAAN_f);
RAAN_rel = 0;

RAAN_dot_init = -((3/2)*(sqrt(mu)*J2*Re^2)/((1-e^2)^2*a^(7/2)))*cos(incl); %Nodal precession in rad/s for initial orbit

posinfo = menu('Propagate spacecraft position (significantly longer simulation time)?', 'Yes', 'No');

%% Iterate through timesteps and track performance over time

%for each time t, use Kepeler's equation to solve for orbital
%position , then convert to perifocal, then to Newtonian
for i = 1:length(t)    
    
    if posinfo == 1
        [posN(:,i)] = PropagateOrbitPos(a, e, incl, anomaly, RAAN, ArgPer, mu, t(i)) ;
    elseif posinfo == 2
        posN(:,i) = zeros(3,1);
    end
    
    %Calculate perturbations in RAAN and ArgPer from J2
    RAAN_dot = -((3/2)*(sqrt(mu)*J2*Re^2)/((1-e^2)^2*a^(7/2)))*cos(incl); %Nodal precession in rad/s
    ArgPer_dot = -((3/2)*(sqrt(mu)*J2*Re^2)/((1-e^2)^2*a^(7/2)))*((5/2)*sin(incl)^2-2); % precession of ArgPer in rad/s
    
    accel = (thrust/sc_mass)/1000; %thrust given in kg*m/s^2, but want accel in km/s^2
    
    if i > 1
        %Calculate new semimajor axis due to acceleration over dt
        %Note: use v and a from previous iteration step
        a0= a;
        a = a/(1-accel*dt(i-1)/v)^2; %Eq 5 from https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-522-space-propulsion-spring-2015/lecture-notes/MIT16_522S15_Lecture6.pdf
%         a = mu/((sqrt(mu/a0)+accel*dt(i-1))^2); %
        alt = a-Re; 
        RAAN = RAAN + RAAN_dot*dt(i-1);
        RAAN_rel = RAAN_rel + (RAAN_dot*dt(i-1)-RAAN_dot_init*dt(i-1));
        ArgPer = ArgPer + ArgPer_dot*dt(i-1);
        delta_v = accel*dt(i-1); 
        v = v + delta_v;
%         a = mu/v^2;


%         [thrust, prev_err, int_err] = pid_altitude(alt,alt_f,alt_0, thrust_max, dt(i-1), int_err, prev_err);
%         [thrust] = sigControl_alt(alt, alt_f,alt_0, thrust_max);
%         [thrust, prev_alt, flag] = AltConroller_simple(alt, alt_f, alt_0, thrust_max, alt_tol, prev_alt, flag);
        [thrust, flag] = RAANController(RAAN_0, RAAN_f, RAAN_rel, thrust_max);

    end
    
    if flag == 1 && resDisp == 0
        fprintf('<strong>Arrived at new orbit!</strong>\n');
        disp(['Transfer time: ', num2str(t(i)/3600), ' hours']);
%         disp(['Delta-V: ' , num2str(v-v_track(1)), ' km/s']); % FOR ALT
        disp(['Delta-V: ' , num2str(sum(abs(delta_v_track(1:i)))), ' km/s']);
        disp(['Propellant expended: ', num2str(m_dot*t(i)), ' kg']);
        disp(['Energy consumed: ', num2str(prop_power*t(i)/3600), ' W-h']);
        disp(['RAAN precession since initiating maneuver: ', num2str(rad2deg(RAAN-RAAN_track(1))), ' deg']);
        resDisp = 1; 
        disp('');
        disp('PLEASE WAIT FOR SIMULATION TO COMPLETE OVER ENTIRE INPUT TIME VECTOR')
    end
    
    RAAN_track(i) = RAAN;
    RAAN_dot_track(i) = RAAN_dot;
    ArgPer_track(i) = ArgPer;
    v_track(i) = v;
    alt_track(i) = a-Re;
    thrust_track(i) = thrust; 
    delta_v_track(i) = delta_v;
    RAAN_rel_track(i) = RAAN_rel;
    
    if flag == 0
        perc_comp = i/length(t)*100;
        disp(['Propagating Orbit... ', num2str(perc_comp), '% Complete']);
    end
    
end

end
