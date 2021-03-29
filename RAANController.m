function [thrust, flag] = RAANController(RAAN_0, RAAN_f, RAAN_rel, thrust_max)
RAAN_diff = RAAN_f-RAAN_0;
RAAN_mid = 0.5*RAAN_diff; 

% If RAAN shift is positive...
if RAAN_diff >= 0
    
    if RAAN_rel <= RAAN_mid %&& RAAN >= RAAN_0 %Remove because orbit will precess to below RAAN_0 before thrust is applied
        thrust = thrust_max;
        flag = 0;
    end

    if RAAN_rel > RAAN_mid && RAAN_rel <= RAAN_diff
        thrust = -thrust_max;
        flag = 0;
    end
    
    if RAAN_rel > RAAN_diff
        thrust = 0; 
        flag = 1;
    end
    
    
end

% If RAAN shift is negative...
if RAAN_diff < 0
    
    if RAAN_rel >= RAAN_mid %&& RAAN_rel <= RAAN
        thrust = -thrust_max;
        flag = 0; 
    end

    if RAAN_rel < RAAN_mid && RAAN_rel >= RAAN_diff
        thrust = thrust_max;
        flag = 0;
    end
    
    if RAAN_rel < RAAN_diff
        thrust = 0; 
        flag = 1;
    end
end
end

% function [thrust, flag] = RAANController(RAAN_0, RAAN_f, RAAN, thrust_max)
% RAAN_diff = RAAN_f-RAAN_0;
% RAAN_mid = RAAN_0 + 0.5*RAAN_diff; 
% 
% % If RAAN shift is positive...
% if RAAN_diff >= 0
%     
%     if RAAN <= RAAN_mid && RAAN >= RAAN_0 %Remove because orbit will precess to below RAAN_0 before thrust is applied
%         thrust = thrust_max;
%         flag = 0;
%     end
% 
%     if RAAN > RAAN_mid && RAAN <= RAAN_f
%         thrust = -thrust_max;
%         flag = 0;
%     end
%     
%     if RAAN > RAAN_f
%         thrust = 0; 
%         flag = 1;
%     end
%     
%     
% end
% 
% % If RAAN shift is negative...
% if RAAN_diff < 0
%     
%     if RAAN >= RAAN_mid && RAAN <= RAAN_0
%         thrust = -thrust_max;
%         flag = 0; 
%     end
% 
%     if RAAN < RAAN_mid && RAAN >= RAAN_f
%         thrust = thrust_max;
%         flag = 0;
%     end
%     
%     if RAAN < RAAN_f
%         thrust = 0; 
%         flag = 1;
%     end
%   
%     
% end
% 
% end
% 

