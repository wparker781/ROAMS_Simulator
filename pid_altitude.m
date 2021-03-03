%% Function Description
%This function can be used as a controller that regulates spacecraft
%attitude by commanding thrust for a non-impulsive maneuver (i.e. low
%thrust). By providing the current altitude, the final desired altitude,
%and the gain constants, we can adjust the commanded thrust such that the
%desired altitude adjustment is effective and does not waste propellant.
%
% NOTE: This function should be used in iteration in combination with
% PropagateOrbit_J2.m. It should be called at every timestep, to determine
% the appropriate thrust command for the following timestep. For best
% results, use a small timestep (<10 minutes). 
%
% NOTE: For simple transfers between two defined orbits in the same plane, this 
% function is not neccesary. 
%
% INPUTS: 
%   alt - The current altitude of the spacecraft [km]
%   alt_f - The desired final altitude [km]
%   kp - proportional gain
%   ki - integral gain
%   kd - derivative gain
%   thrust_max - maximum thrust from propulsion system
%   dt - timestep for simulation per iteration in PropagateOrbit_J2 [s]
%   int_err - tracks an approximation of the integrated error
%   prev_err - the error recorded in the previous timestep (for
%   approximating the derivative of error). 
%   
% OUTPUTS:
%   thrust - thrust command for current timestep
%   err - the error between current and desired altitude [km]
%   int_err - the approximation of the integrated error
%
% William Parker - March 2021
%% 
function [thrust, err, int_err] = pid_altitude(alt,alt_f,alt_0, thrust_max, dt, int_err, prev_err)
    kp = 1;
    ki = 1; 
    kd = 1;

    err = alt_f-alt;
    norm_err = err/abs(alt_f-alt_0);% guarantees that error will be [-1,1] (won't exceed max thrust).
    int_err = int_err + norm_err*dt;
    % deriv_err = (prev_err-err)/dt;
    thrust = thrust_max*(kp*norm_err); %+ ki*(int_err)+ kd*deriv_err);

end



