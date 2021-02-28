% Repeating Ground Track Calculator
% Based on: http://www.aiaahouston.org/Horizons/RepeatGTr0.pdf
% Written by Daniel Miller, 10/6/2020, for 16.851
% Apologies to anyone using this, it's very quick and dirty

clear
clc
close all

% ----------------------------------------------------------------------- %
% Set the desired Apogee, Perigee, inclination
set = 566.897;

Ha = set;       % Desired Apogee Altitude [km]
Hp = set;       % Desired Perigee Altitude [km]
Hfix = Hp;      % Required Altitude [km] -- may by Hp or Ha, but it hit it
inc = 100.965;       % Required inclination
% ----------------------------------------------------------------------- %

% Orbital Parameters
mu = 398600.4415;
R = 6378.1363;
J2 = 0.0010826269;
w = 7.292115e-05;

% magic (orbital mechanics stuff)
a = R + 0.5*(Ha + Hp);
e = 1 - (R+Hp)/a;
n = 1/sqrt(a^3/ mu);
Pk = 2*pi*sqrt(a^3/mu);

Pomega = Pk*(1-1.5*J2*(R/a)^2 * (3 - 4*sind(inc)^2));

p = a*(1-e^2);

omegaDot = -1.5*(n*R^2 * J2)/(p^2) * cosd(inc);
lambdaDot = omegaDot - w;
k = linspace(1,20,20);
sigma = (k.*Pomega.*lambdaDot)./(2*pi);

% sigma needs to be an integer. Adjust semi-major axis using fsolve
remainders = abs(sigma - round(sigma));
[sortedRemainders,I] = sort(remainders);

params.mu = mu;
params.J2 = J2;
params.R = R;
params.inc = inc;
params.w = w;
params.n = n;
params.e = e;

solutions = zeros(4,10);
% solverOpts = optimset('TolFun',1e-08,'TolX',1e-10);
solverOpts = optimoptions('fsolve');
solverOpts.OptimalityTolerance = 1e-10;
solverOpts.StepTolerance = 1e-10;
for j = 1:10
    initGuess = Ha;
    [x,fval, exitflag, ~] = fsolve(@(Hfree)fixSigma(Hfree, Hfix, k(I(j)), params),initGuess,solverOpts);
    if exitflag > 0
        solutions(1,j) = x;
        solutions(2,j) = Hfix;
        solutions(3,j) = abs(((x + R) - (Hfix + R))/((x + R) + (Hfix + R)));
        solutions(4,j) =  k(I(j));
    else
        solutions(3,j) = 100;
    end
end

% Pick the most circular orbit from list of solutions
[tab, solution] = generateTable(solutions,params);
disp(tab)

% Propogate orbit using J2 model, convert to lat/long/alt, account for the
% rotation of the Earth, and finally plot the results
a = (solution.Ha + solution.Hp)/2 + R;
vInit = sqrt(mu*(2/(Hfix + R) - 1/a));

% initState = [Hfix + R; 0; 0; 0; vInit*cosd(inc); vInit*sind(inc)];
% options = odeset('RelTol',1e-12,'AbsTol',1e-12);
% [t,y] = ode45(@(t,x)j2Dynamics(t,x,a),[0,solution.T*60],initState,options);

% vInit = tweakVel(solution,Hfix, params);

wgs84 = wgs84Ellipsoid('kilometers');
initState = [Hfix + R; 0; 0; 0; vInit*cosd(inc); vInit*sind(inc)];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,y] = ode45(@(t,x)j2Dynamics(t,x,a),[0,86400*1],initState,options);

lla = ecef2lla(y(:,1:3).*1000);
lla(:,2) = lla(:,2) - (w.*t)*(180/pi);
plotGroundTracks(lla)

% tweak(solution, Hfix, params)


% ----------------------------------------------------------------------- %
function sigmaRemainder = fixSigma(Hfree, Hfix, k, params)
%fixSigma Summary of this function goes here
%   Detailed explanation goes here

    mu = params.mu;
    J2 = params.J2;
    R = params.R;
    inc = params.inc;
    w = params.w;
%     n = params.n;
    
    a = R + 0.5*(Hfree + Hfix);
    if Hfree < Hfix
        e = 1 - (R+Hfree)/a;
    else
        e = 1 - (R+Hfix)/a;
    end
    n = 1/sqrt(a^3/ mu);
    Pk = 2*pi*sqrt(a^3/mu);
    Pomega = Pk*(1-1.5*J2*(R/a)^2 * (3 - 4*sind(inc)^2));
    p = a*(1-e^2);
    omegaDot = -1.5*(n*R^2 * J2)/(p^2) * cosd(inc);
    lambdaDot = omegaDot - w;
    sigma = (k.*Pomega.*lambdaDot)./(2*pi);
    sigmaRemainder = abs(sigma - round(sigma));
end

function [dX] = j2Dynamics(~,X,a)
% j2Dynamics function to be integrated to produce J2 orbits
    J2 = 0.0010826269;
    mu = 398600.4415;

    x = X(1);
    y = X(2);
    z = X(3);
    dx = X(4);
    dy = X(5);
    dz = X(6);

    r = norm(X(1:3));

    dX(1) = dx;
    dX(2) = dy;
    dX(3) = dz;
    dX(4) = -(mu*x/r^3)*(1 + (3/2)*J2*(a^2/r^2) - (15/2)*J2*(a^2*z^2)/r^4);
    dX(5) = -(mu*y/r^3)*(1 + (3/2)*J2*(a^2/r^2) - (15/2)*J2*(a^2*z^2)/r^4);
    dX(6) = -(mu*z/r^3)*(1 + (9/2)*J2*(a^2/r^2) - (15/2)*J2*(a^2*z^2)/r^4);
    dX = dX';
end

function plotGroundTracks(pos)
% plotGroundTracks
    lat = pos(:,1);
    lon = wrapTo180(pos(:,2));
    mstruct = defaultm('mercator');
    load coastlines
    mstruct.geoid = referenceEllipsoid('wgs84','kilometers');
%     mstruct.maplatlimit = [-80, 80];
    mstruct.maplatlimit = [-90, 90];

    mstruct.maplonlimit = [-180, 180];
    [scLat, scLong] = maptriml(lat,lon, ...
     mstruct.maplatlimit,mstruct.maplonlimit);
%     [x,y] = projfwd(mstruct,scLat,scLong);
    figure
%     fig = plot(x,y);
    geoshow(scLat,scLong)
   
    hold on
    [latt,lont] = maptriml(coastlat,coastlon, ...
     mstruct.maplatlimit,mstruct.maplonlimit);
%     [x2,y2] = projfwd(mstruct,latt,lont);
%     plot(x2,y2)
    geoshow(latt,lont,'Color','k')
    
%     title('Ground Track','Fontsize',14)
    xlabel('Longitude [deg]','Fontsize',14)
    ylabel('Latitude [deg]','Fontsize',14)
    legend('Ground Track','Fontsize',14)
    grid on
%     axis equal
    
end

function [tab, solution] = generateTable(solutions,params)
    [sorted_e,ind] = sort(solutions(3,:));
    sortedSolns = solutions(:,ind);
    good_e = find(sorted_e ~= 100);
    good_hp = find(sortedSolns(1,:) > 0);
    good_ha = find(sortedSolns(2,:) > 0);
    goodIndices = intersect(good_e,intersect(good_hp,good_ha));
    
    varNames = {'Perigee [km]','Apogee [km]','Eccentricity',...
            'Inclination [deg]', 'Period [min]','Orbit Repeat #',...
            'Sigma Remainder'};
    tableArray = zeros(length(goodIndices),7);
    
    inc = params.inc;
    R = params.R;
    mu = params.mu;
    
    for j = 1:length(goodIndices)
        i = goodIndices(j);
        T = 2*pi*sqrt((mean(sortedSolns(1:2,i))+R)^3 / mu)/60;
        sigmaRemainder = fixSigma(sortedSolns(1,i), sortedSolns(2,i),...
            sortedSolns(4,i), params);
        if sortedSolns(1,i) < sortedSolns(2,i)
            tableArray(j,:) = [sortedSolns(1,i),sortedSolns(2,i),...
                sortedSolns(3,i),inc,T,sortedSolns(4,i),sigmaRemainder];
        else
            tableArray(j,:) = [sortedSolns(2,i),sortedSolns(1,i),...
                sortedSolns(3,i),inc,T,sortedSolns(4,i),sigmaRemainder];
        end
    end
    
    solution.Hp = tableArray(1,1);
    solution.Ha = tableArray(1,2);
    solution.e = tableArray(1,3);
    solution.inc = tableArray(1,4);
    solution.T = tableArray(1,5);
    solution.k = tableArray(1,6);
    tab = array2table(tableArray,'VariableNames',varNames); 
end

% function vInitP = tweakVel(initSoln, Hfix, params)
% % hair brained scheme that didn't really work
%     a = (initSoln.Ha + initSoln.Hp)/2 + params.R;
%     vInit = sqrt(params.mu*(2/(Hfix + params.R) - 1/a));
%     initState = [Hfix + params.R; 0; 0; 0; vInit*cosd(params.inc); vInit*sind(params.inc)];
%     options = odeset('RelTol',1e-12,'AbsTol',1e-12);
%     [~,y] = ode45(@(t,x)j2Dynamics(t,x,a),[0,initSoln.T*60],initState,options);
%     alt = vecnorm(y(:,1:3),2,2) - params.R;
%     
%     if Hfix == initSoln.Ha
%         % then lets fix Hp
%         Hp_true = min(alt);
%         if Hp_true < initSoln.Hp
%             for perturb = 1:0.0001:1.01
%                 vInitP = vInit*perturb;
%                 initState = [initSoln.Ha + params.R; 0; 0; 0;...
%                     vInitP*perturb*cosd(params.inc); vInitP*perturb*sind(params.inc)];
%                 [~,y] = ode45(@(t,x)j2Dynamics(t,x,a),[0,initSoln.T*60],initState,options);
%                 alt = vecnorm(y(:,1:3),2,2) - params.R;
%                 if abs(initSoln.Hp - min(alt))/initSoln.Hp < 0.005
%                     break 
%                 end
%             end
%         else
%             for perturb = linspace(1,0.99,101)
%                 vInitP = vInit*perturb;
%                 initState = [initSoln.Ha + params.R; 0; 0; 0;...
%                     vInitP*perturb*cosd(params.inc); vInitP*perturb*sind(params.inc)];
%                 [~,y] = ode45(@(t,x)j2Dynamics(t,x,a),[0,initSoln.T*60],initState,options);
%                 alt = vecnorm(y(:,1:3),2,2) - params.R;
%                 if abs(initSoln.Hp - min(alt))/initSoln.Hp < 0.005
%                     break 
%                 end
%             end
%         end
%     else
%         % then lets fix Ha
%         Ha_true = max(alt);
%         if Ha_true < initSoln.Ha
%             for perturb = 1:0.0001:1.01
%                 vInitP = vInit*perturb;
%                 initState = [initSoln.Hp + params.R; 0; 0; 0;...
%                     vInitP*perturb*cosd(params.inc); vInitP*perturb*sind(params.inc)];
%                 [~,y] = ode45(@(t,x)j2Dynamics(t,x,a),[0,initSoln.T*60],initState,options);
%                 alt = vecnorm(y(:,1:3),2,2) - params.R;
%                 disp(abs(initSoln.Ha - max(alt))/initSoln.Ha)
%                 if abs(initSoln.Ha - max(alt))/initSoln.Ha < 0.005
%                     break 
%                 end
%             end
%         else
%             for perturb = linspace(1,0.99,101)
%                 vInitP = vInit*perturb;
%                 initState = [initSoln.Hp + params.R; 0; 0; 0;...
%                     vInitP*perturb*cosd(params.inc); vInitP*perturb*sind(params.inc)];
%                 [~,y] = ode45(@(t,x)j2Dynamics(t,x,a),[0,initSoln.T*60],initState,options);
%                 alt = vecnorm(y(:,1:3),2,2) - params.R;
%                 if abs(initSoln.Ha - max(alt))/initSoln.Ha < 0.005
%                     break 
%                 end
%             end
%         end
%     end
%     
% end
% 
% function tweak(initSoln, fixed, params)
% % tweak. Don't trust fsolve? This will "poke" the solution by +/- 5% to
% %   check the solution, just in case.
%     if initSoln.Ha == fixed
%         notFixed_ = initSoln.Hp;
%     else
%         notFixed_ = initSoln.Ha;
%     end
%     
%     noise = 0.95:0.001:1.05;
%     data = zeros(101,2);
% 
%     for i = 1:length(noise)
%         notFixed = notFixed_*noise(i);
%         sigmaRemainder = fixSigma(notFixed, fixed, initSoln.k, params);
%         data(i,1) = notFixed;
%         data(i,2) = sigmaRemainder;
%     end
%     
%     [M,I] = min(data(:,2));
%     disp(data(I,1))
%     disp(M)
%     figure
%     plot(data(:,1),data(:,2))
% end
