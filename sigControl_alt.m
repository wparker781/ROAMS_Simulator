%% Function Description
% Uses a sigmoid to approximate a time-optimal low thrust maneuver
% (bang-bang pos max to low max) with spacecraft rotation time for
% repositioning thrusters (not necessary for simple transfers, which only
% apply thrust in one direction). 
%
% William Parker - March 2021
%%
function [thrust] = sigControl_alt(alt, alt_f,alt_0, max_thrust)
err = alt_f-alt;
norm_err = err/abs(alt_f-alt_0);% guarantees that error will be [-1,1] (won't exceed max thrust).
%define sigmoid function that will be used to approximate optimal trajectory
sig_const = 30; %a value of 30 makes the sigmoid slope very large w/ sharp changes (which is what we're looking for). 
if norm_err > 0 % if we're trying to increase altitude...
    if norm_err-0.5 < tol % avoid getting stuck with no control input
        s = 2*(1./(1+exp(-sig_const.*(norm_err-0.5)))) - 1;
    
else
    s = 2*(1./(1+exp(-sig_const.*(norm_err+0.5)))) - 1;
end

thrust = s*max_thrust; 

end




