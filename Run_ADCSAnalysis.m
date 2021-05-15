%% Function Description
% input thrust, thrust angle, COM offset. Calculate momentum storage
% required, frequency of desaturation required, whether this is within spec
% for the reaction wheel assembly and magnetorquer from AssignParams.m. 
% 
% Working principle: T = Ia. Torque from off-center thrusting should be as
% small as possible. Any torque produced from this misalignment will need
% to be corrected for by the reaction wheel assembly (RWA). The RWA can
% oppose the torque from the thruster misalignment by accelerating the
% reaciton wheels, but they can only accelerate to a certain maximum speed.
% Beyond that wheel speed, the reation wheels need to be desaturated with 
% magnetorquers to retain control of the spacecraft. 

%% ASSIGN PARAMETERS
thrust_level = 2; 
[planet,t,orbInit,orbFin,numHT,numLT,sc,prop, adcs] = AssignParams(thrust_level);

save('adcs.mat', '-struct', 'adcs');
save('sc.mat', '-struct', 'sc');
save('prop.mat', '-struct', 'prop');

load('adcs.mat');
load('sc.mat');
load('prop.mat');

%% Perform analysis for a variety of offset distances
offset = [0:0.0001:0.05];
T_thrust = thrust.*offset;

% assume magnetic field strength for the desired orbit varies between
% 20,000-45000 nT.
% (https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml?useFullSite=true#igrfwmm)
% Take average for B
B = ((20000+45000)/2)/10^9;
% Assume NCTR-M002 Magnetorquer Rod (cubesat scale, but not necessarily
% what is used in BCT bus)
mag_mom = 0.2; %Am^2
% assume ideal geometry (B and mag_mom orthogonal)
T_mag = B*mag_mom; %torque in N-m. 

%to counter-balance T_thrust, need to match with the reaction wheel
%assembly
T_rwa = zeros(length(offset),1);
t_desat = zeros(length(offset),1); 
for i =  1: length(offset)
    if T_thrust(i) < maxT
        T_rwa(i) = T_thrust(i);
    else 
        T_rwa(i) = 0;
    end
    t_desat(i) = momStore/T_rwa(i); %time between desaturation in s
end
   
%% Create plots showing time between desaturation for different offset distances 

semilogy(offset*100, t_desat/3600, 'r-')
xlabel('Offset distance [cm]');
ylabel('Time until desaturation necessary [hours]');
title('Desaturation Interval assuming RWP050 RWA')
hold on

% Now, calculate necessary RWA torque assuming that it's getting some help
% from magnetorquers through the duration of the burn
T_rwa = zeros(length(offset),1);
t_desat = zeros(length(offset),1); 
for i =  1: length(offset)
    if T_thrust(i) < maxT + T_mag
        T_rwa(i) = T_thrust(i)-T_mag;
    else 
        T_rwa(i) = 0;
    end
    if T_rwa < 0
        T_rwa = 0;
    end
    t_desat(i) = momStore/T_rwa(i); %time between desaturation in s
end

semilogy(offset*100, t_desat/3600, 'r--')
% legend('Enpulsion Nano (Alone)', 'Enpulsion Nano (with Mag T)');
xlim([0,max(offset)*100]);

hold on
