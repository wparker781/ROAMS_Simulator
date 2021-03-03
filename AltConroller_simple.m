%% Function Description
% Simple controller for commanding thrust for a low thrust transfer where
% the initial and final orbit are well defined and lie in the same plane. 
%
% William Parker - March 2021
%%
function [thrust, prev_alt, flag] = AltConroller_simple(alt, alt_f, alt_0, thrust_max, alt_tol, prev_alt, flag)

if flag == 1
    thrust = 0; 
    prev_alt = alt;
    flag = 1;
    return
end

if (alt-alt_f) > alt_tol
    thrust = -thrust_max;
else if (alt-alt_f) < -alt_tol
        thrust = thrust_max;
    else if abs(alt-alt_f) <= alt_tol
            thrust = 0; 
            
        end
    end
end

%If you've crossed over the desired altitude, stop applying thrust (switch
%over to fine maneuvering). 
if (alt-alt_f)*(prev_alt-alt_f) < 0
    thrust = 0; 
    %throw a flag to indicate that we don't need to apply more thrust now
    flag = 1; 
end

prev_alt = alt; 

end
