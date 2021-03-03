function [posN,RAAN_track,RAAN_dot_track,ArgPer_track, v_track, alt_track, thrust_track, t] = PropagateOrbit_J2(t,mu,a,e,h,incl,RAAN,ArgPer,anomaly, thrust_max, sc_mass, m_dot,prop_power, r_f, alt_tol, Re)
%% Function Description
% This function is responsible for taking orbit parameters (including actual 
% anomaly) as inputs, and outputting a position vector at each time interval t 
% from the center of the attracting body to the object of interest based on the
% orbit parameters. 
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

%% Initialize variables
syms E theta
solE = zeros(length(t),1);
theta = zeros(length(t),1);
r = zeros(length(t),1);
PosOrb = zeros(3,length(t));
posPF = zeros(3,length(t));
RAAN_track = zeros(1,length(t));
RAAN_dot_track = zeros(1,length(t)); 
ArgPer_track = zeros(1,length(t));
v_track = zeros(1,length(t)); 
alt_track = zeros(1,length(t));
thrust_track = zeros(1,length(t)); 
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

%% Iterate through timesteps and track performance over time

%for each time t, use Kepeler's equation to solve for orbital
%position , then convert to perifocal, then to Newtonian
for i = 1:length(t)    
    %Create direction cosine matrix for converting a newtonian frame to the
    %perifocal
    C_N2PF = [-sin(RAAN)*cos(incl)*sin(ArgPer)+cos(RAAN)*cos(ArgPer), cos(RAAN)*cos(incl)*sin(ArgPer)+sin(RAAN)*cos(ArgPer), sin(incl)*sin(ArgPer);...
                -sin(RAAN)*cos(incl)*cos(ArgPer)-cos(RAAN)*sin(ArgPer), cos(RAAN)*cos(incl)*cos(ArgPer)-sin(RAAN)*sin(ArgPer), sin(incl)*cos(ArgPer);...
                sin(RAAN)*sin(incl), -cos(RAAN)*sin(incl), cos(incl)];
    % Kep(i) = sqrt(mu/a)*(t(i)) == e*sin(E)-E;
    solE(i) = vpasolve(sqrt(mu/a^3)*(t(i)) == E-e*sin(E),E);
    theta(i) = 2*atan(sqrt((1+e)/(1-e))*tan(solE(i)/2))+anomaly;
    % solTheta(i) = vpasolve(tan(theta/2)==sqrt((1+e)/(1-e))*tan(solE(i)/2),theta);
    r(i) = h^2/mu*(1+e*cos(theta(i))).^-1;
    PosOrb(:,i) = [r(i);0;0];
    C_orb2PF = [cos(theta(i)) sin(theta(i)) 0;-sin(theta(i)) cos(theta(i)) 0;0 0 1];
    posPF(:,i) = ((C_orb2PF)'*PosOrb(:,i));
    posN(:,i) = C_N2PF'*posPF(:,i);
    
    %Calculate perturbations in RAAN and ArgPer from J2
    RAAN_dot = -((3/2)*(sqrt(mu)*J2*Re^2)/((1-e^2)^2*a^(7/2)))*cos(incl); %Nodal precession in rad/s
    ArgPer_dot = -((3/2)*(sqrt(mu)*J2*Re^2)/((1-e^2)^2*a^(7/2)))*((5/2)*sin(incl)^2-2); % precession of ArgPer in rad/s
    
    accel = (thrust/sc_mass)/1000; %thrust given in kg*m/s^2, but want accel in km/s^2
    
    if i > 1
        %Calculate new semimajor axis due to acceleration over dt
        %Note: use v and a from previous iteration step
        a = a/(1-accel*dt(i-1)/v)^2; %Eq 5 from https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-522-space-propulsion-spring-2015/lecture-notes/MIT16_522S15_Lecture6.pdf
        alt = a-Re; 
        RAAN = RAAN + RAAN_dot*dt(i-1);
        ArgPer = ArgPer + ArgPer_dot*dt(i-1);
        delta_v = accel*dt(i-1); 
        v = v + delta_v;
%         a = mu/v^2;


%         [thrust, prev_err, int_err] = pid_altitude(alt,alt_f,alt_0, thrust_max, dt(i-1), int_err, prev_err);
%         [thrust] = sigControl_alt(alt, alt_f,alt_0, thrust_max);
        [thrust, prev_alt, flag] = AltConroller_simple(alt, alt_f, alt_0, thrust_max, alt_tol, prev_alt, flag);
    
    end
    
    if flag == 1 && resDisp == 0
        fprintf('<strong>Arrived at new orbit!</strong>\n');
        disp(['Transfer time: ', num2str(t(i)/3600), ' hours']);
        disp(['Delta-V: ' , num2str(v-v_track(1)), ' km/s']);
        disp(['Propellant expended: ', num2str(m_dot*t(i)), ' kg']);
        disp(['Energy consumed: ', num2str(prop_power*t(i)/3600), ' W-h']);
        disp(['RAAN precession since initiating maneuver: ', num2str(rad2deg(RAAN-RAAN_track(1))), ' deg']);
        resDisp = 1; 
    end
    
    RAAN_track(i) = RAAN;
    RAAN_dot_track(i) = RAAN_dot;
    ArgPer_track(i) = ArgPer;
    v_track(i) = v;
    alt_track(i) = a-6378.15;
    thrust_track(i) = thrust; 
    
    if flag == 0
        perc_comp = i/length(t)*100;
        disp(['Propagating Orbit... ', num2str(perc_comp), '% Complete']);
    end
    
end

end
