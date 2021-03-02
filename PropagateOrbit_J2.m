function [posN,RAAN_track,RAAN_dot_track, v_track, a_track, t] = PropagateOrbit_J2(t,mu,a,e,h,incl,RAAN,ArgPer,anomaly, thrust, sc_mass)
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
% maneuvering. For this orbit propagation tool, we assume continuous thrust
% over the duration of the input time vector (N). 

% OUTPUTS:
% PosN = an nx3 matrix of position vectors with 3 components for each time t [km]
% RAAN (updated at each timestep because of the J2 perturbation)
% t = nx1 array of times since initial actual anomaly (s)

% NOTE: Assuming Earth for radius and J2
% 
% Developed by William Parker, November 2019
% 3/1/2021 - Modified to add J2 perturbation 

%% Initialize variables
syms E theta
solE = zeros(length(t),1);
theta = zeros(length(t),1);
r = zeros(length(t),1);
PosOrb = zeros(3,length(t));
posPF = zeros(3,length(t));
RAAN_track = zeros(1,length(t));
RAAN_dot_track = zeros(1,length(t)); 
v_track = zeros(1,length(t)); 

J2 = 1.08263e-3; 
dt = diff(t);

%convert angles in degrees to radians
incl = deg2rad(incl);
RAAN = deg2rad(RAAN);
ArgPer = deg2rad(ArgPer);
anomaly = deg2rad(anomaly);

%Calculate velocity and accel of sc in the initial orbit
v = sqrt(mu/a);
accel = (thrust/sc_mass)/1000; %thrust given in kg*m/s^2, but want accel in km/s^2


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
    RAAN_dot = -((3/2)*(sqrt(mu)*J2*6378.15^2)/((1-e^2)^2*a^(7/2)))*cos(incl); %Nodal precession in rad/s
    ArgPer_dot = -((3/2)*(sqrt(mu)*J2*6378.15^2)/((1-e^2)^2*a^(7/2)))*((5/2)*sin(incl)^2-2); % precession of ArgPer in rad/s
    
    %Calculate change in semi-major axis due to low-thrust during timestep
    %Assume that for the duration of the burn, we retain a circular orbit
    %with 0 eccentricity. 
%     accel = accel + thrust/sc_mass; 
    
   
    if i > 1
        %Calculate new semimajor axis due to acceleration over dt
        %Note: use v and a from previous iteration step
        a = a/(1-accel*dt(i-1)/v)^2; %Eq 5 from https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-522-space-propulsion-spring-2015/lecture-notes/MIT16_522S15_Lecture6.pdf
        RAAN = RAAN + RAAN_dot*dt(i-1);
        ArgPer = ArgPer + ArgPer_dot*dt(i-1);
        delta_v = accel*dt(i-1); 
        v = v + delta_v;
%         a = mu/v^2;
       
    end
    
    RAAN_track(i) = RAAN;
    RAAN_dot_track(i) = RAAN_dot; 
    v_track(i) = v;
    a_track(i) = a;
    
            
    perc_comp = i/length(t)*100;
    disp(['Propagating Orbit... ', num2str(perc_comp), '% Complete']);
    
end

end