%% Function Description
% Input Keplerian orbit elements (without true anomaly), output the orbital position in the
% Earth-centered newtonian frame at the future time provided (in s).
%
% William Parker - March 2021
% 
function [posN] = PropagateOrbitPos(a, e, incl,anomaly, RAAN, ArgPer, mu, t)    

    
    %Create direction cosine matrix for converting a newtonian frame to the
    %perifocal
    h = sqrt(a*mu);
    syms E theta
    C_N2PF = [-sin(RAAN)*cos(incl)*sin(ArgPer)+cos(RAAN)*cos(ArgPer), cos(RAAN)*cos(incl)*sin(ArgPer)+sin(RAAN)*cos(ArgPer), sin(incl)*sin(ArgPer);...
                -sin(RAAN)*cos(incl)*cos(ArgPer)-cos(RAAN)*sin(ArgPer), cos(RAAN)*cos(incl)*cos(ArgPer)-sin(RAAN)*sin(ArgPer), sin(incl)*cos(ArgPer);...
                sin(RAAN)*sin(incl), -cos(RAAN)*sin(incl), cos(incl)];
    % Kep(i) = sqrt(mu/a)*(t(i)) == e*sin(E)-E;
    if t > 0
        solE = vpasolve(sqrt(mu/a^3)*(t) == E-e*sin(E),E);
    elseif t == 0
        solE = sqrt(mu/a^3)*t; %assumes circular orbits!!
    else
        disp('t must be >= 0');
    end
    
    theta = 2*atan(sqrt((1+e)/(1-e))*tan(solE/2))+anomaly;
    % solTheta(i) = vpasolve(tan(theta/2)==sqrt((1+e)/(1-e))*tan(solE(i)/2),theta);
    r = h^2/mu*(1+e*cos(theta)).^-1;
    PosOrb = [r;0;0];
    C_orb2PF = [cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1];
    posPF = ((C_orb2PF)'*PosOrb);
    posN = C_N2PF'*posPF;
    
end

%% For reference, earlier iterative implementation 

%     %Create direction cosine matrix for converting a newtonian frame to the
%     %perifocal
%     C_N2PF = [-sin(RAAN)*cos(incl)*sin(ArgPer)+cos(RAAN)*cos(ArgPer), cos(RAAN)*cos(incl)*sin(ArgPer)+sin(RAAN)*cos(ArgPer), sin(incl)*sin(ArgPer);...
%                 -sin(RAAN)*cos(incl)*cos(ArgPer)-cos(RAAN)*sin(ArgPer), cos(RAAN)*cos(incl)*cos(ArgPer)-sin(RAAN)*sin(ArgPer), sin(incl)*cos(ArgPer);...
%                 sin(RAAN)*sin(incl), -cos(RAAN)*sin(incl), cos(incl)];
%     % Kep(i) = sqrt(mu/a)*(t(i)) == e*sin(E)-E;
%     solE(i) = vpasolve(sqrt(mu/a^3)*(t(i)) == E-e*sin(E),E);
%     theta(i) = 2*atan(sqrt((1+e)/(1-e))*tan(solE(i)/2))+anomaly;
%     % solTheta(i) = vpasolve(tan(theta/2)==sqrt((1+e)/(1-e))*tan(solE(i)/2),theta);
%     r(i) = h^2/mu*(1+e*cos(theta(i))).^-1;
%     PosOrb(:,i) = [r(i);0;0];
%     C_orb2PF = [cos(theta(i)) sin(theta(i)) 0;-sin(theta(i)) cos(theta(i)) 0;0 0 1];
%     posPF(:,i) = ((C_orb2PF)'*PosOrb(:,i));
%     posN(:,i) = C_N2PF'*posPF(:,i);

