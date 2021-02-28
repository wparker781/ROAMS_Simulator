% Repeating Ground Track Calculator
% Based on: http://www.aiaahouston.org/Horizons/RepeatGTr0.pdf
% Written by Daniel Miller, 10/6/2020, for 16.851
% Apologies to anyone using this, it's very quick and dirty

%% RGT Calculation
clear
clc
close all

% ----------------------------------------------------------------------- %
% Set the desired Apogee, Perigee, inclination
alt = 400;
Hfix = alt;      % Required Altitude [km] -- may by Hp or Ha, but it hit it
inc = 51.6;     % Required inclination
% ----------------------------------------------------------------------- %

% Orbital Parameters
mu = 398600.4415;
R = 6378.1363;
J2 = 0.0010826269;
w = 7.292115e-05;

% magic (orbital mechanics stuff)
a = R + alt;
e = 0;
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

solutions = zeros(3,10);
solverOpts = optimoptions('fsolve');
solverOpts.OptimalityTolerance = 1e-10;
solverOpts.StepTolerance = 1e-10;
for j = 1:10
    [x,fval, exitflag, ~] = fsolve(@(alt_)fixSigma(alt_, k(I(j)), params),alt,solverOpts);
    if exitflag > 0
        solutions(1,j) = x;
        solutions(2,j) = fval;
%         solutions(2,j) = abs(((x + R) - (Hfix + R))/((x + R) + (Hfix + R)));
        solutions(3,j) =  k(I(j));
    else
        solutions(3,j) = 100;
    end
end

% Pick the most circular orbit from list of solutions
[tab, solution] = generateTable(solutions,params);
disp(tab)

% Propogate orbit using J2 model, convert to lat/long/alt, account for the
% rotation of the Earth, and finally plot the results
vInit = sqrt(mu/solution.a);

wgs84 = wgs84Ellipsoid('kilometers');
initState = [solution.a; 0; 0; 0; vInit*cosd(inc); vInit*sind(inc)];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,y] = ode45(@(t,x)j2Dynamics(t,x,a),[0:86400*1],initState,options);

lla = ecef2lla(y(:,1:3).*1000);
lla(:,2) = lla(:,2) - (w.*t)*(180/pi);

plotGroundTracks(lla,'other')

% vInitGOM = sqrt(mu/(547.5+R));
% initStateGOM = [(547.5+R); 0; 0; 0; vInitGOM*cosd(inc); vInitGOM*sind(inc)];

% vInitGOM = sqrt(mu/(400+R));
% initStateGOM = [(400+R); 0; 0; 0; vInitGOM*cosd(inc); vInitGOM*sind(inc)];
% 
% [tGOM,yGOM] = ode45(@(t,x)j2Dynamics(t,x,a),[0,86400*7],initStateGOM,options);
% llaGOM = ecef2lla(yGOM(:,1:3).*1000);
% llaGOM(:,2) = llaGOM(:,2) - (w.*tGOM)*(180/pi);
% plotGroundTracks(llaGOM,'other')


%% Calculate RAAN for a given ground target -------------------------------
clc
close all

targetLat = 42.3601;
targetLong = -71.0942;

fprintf('The target''s coordinates are %3.4f,%3.4f\n\n',targetLat,targetLong)

% calculate inite specific angular momentum, RAAN
h = cross(y(1,1:3),y(1,4:6));
n = cross([0,0,1],h);
if n(2) >= 0
    omega = acosd(n(1)/norm(n));
else
    omega = 360 - acosd(n(1)/norm(n));
end
fprintf('The original RAAN is %3.4f\n\n',omega)

% find points with the correct latitude
[~,LatIndices] = sort(abs(lla(:,1)-targetLat));

% find points cloest in longitude
[~,LongIndicesOfLat] = sort(abs(lla(LatIndices(1:50),2)-targetLong));
bestIndex = LatIndices(LongIndicesOfLat(1));
fprintf('The initial closest point near that latitude is %3.4f, %3.4f\n\n',...
   lla(bestIndex,1),lla(bestIndex,2));

deltaOmega =(targetLong - lla(bestIndex,2));
reqdOmega = omega + deltaOmega;
fprintf('A viable RAAN is %3.4f\n\n',reqdOmega)

R = [cosd(deltaOmega),  sind(deltaOmega),   0;
     -sind(deltaOmega), cosd(deltaOmega),   0;
     0,                 0,                  1];

initPos= initState(1:3)' * R;
initVel= initState(4:6)' * R;
initStateNew = [initPos, initVel]';

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[tNew,yNew] = ode45(@(t,x)j2Dynamics(t,x,a),[0,86400*1],initStateNew,options);
llaNew = ecef2lla(yNew(:,1:3).*1000);
w = 7.292115e-05;
llaNew(:,2) = llaNew(:,2) - (w.*tNew)*(180/pi);

% calculate inite specific angular momentum, RAAN
hNew = cross(yNew(1,1:3),yNew(1,4:6));
nNew = cross([0,0,1],hNew);
if nNew(2) >= 0
    omegaNew = acosd(nNew(1)/norm(nNew));
else
    omegaNew = 360 - acosd(nNew(1)/norm(nNew));
end
fprintf('The new RAAN is %3.4f\n\n',omegaNew)


% find points with the correct latitude
[~,LatIndicesNew] = sort(abs(llaNew(:,1)-targetLat));

% find points cloest in longitude
[~,LongIndicesOfLatNew] = sort(abs(llaNew(LatIndicesNew(1:50),2)-targetLong));
bestIndexNew = LatIndicesNew(LongIndicesOfLatNew(1));
fprintf('The new closest point near that latitude is %3.4f, %3.4f\n\n',...
   llaNew(bestIndexNew,1),llaNew(bestIndexNew,2));

% eVect = cross(yNew(1,4:6),cross(yNew(1,1:3),yNew(1,4:6)))./mu - yNew(1,1:3)./norm(yNew(1,1:3));
% argPer = acosd(dot(nNew,eVect)/norm(nNew)*norm(eVect));
% 
% fprintf('The argument of periapsis is %3.4f\n\n',argPer)

%% GOM overfly visit frequency
ecefGOM = lla2ecef(llaGOM);
ecefTarget = lla2ecef([targetLat,targetLong,0]);
dist2target = vecnorm(ecefGOM - ecefTarget,2,2);
figure
hold on
plot(tGOM,dist2target)
plot(tGOM,ones(length(tGOM),1).*739.3749*1000);

%% Functions ----------------------------------------------------------- %
function sigmaRemainder = fixSigma(alt, k, params)
%fixSigma Summary of this function goes here
%   Detailed explanation goes here

    mu = params.mu;
    J2 = params.J2;
    R = params.R;
    inc = params.inc;
    w = params.w;

    e = 0;
    a = R + alt;
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

function plotGroundTracks(pos,orbit)
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
    
    if strcmp(orbit,'ROM')
        figure
        geoshow(scLat,scLong,'Color','red','Linewidth',2)

        hold on
        [latt,lont] = maptriml(coastlat,coastlon, ...
         mstruct.maplatlimit,mstruct.maplonlimit);
        geoshow(latt,lont,'Color','k')

        xlabel('Longitude [deg]','Fontsize',14)
        ylabel('Latitude [deg]','Fontsize',14)
        legend('Ground Track','Fontsize',14)
        grid on
    elseif strcmp(orbit,'GOM')
        geoshow(scLat,scLong)
        legend('RGT Ground Track','GOM Ground Track','Fontsize',14)
    else
        figure
        [latt,lont] = maptriml(coastlat,coastlon, ...
         mstruct.maplatlimit,mstruct.maplonlimit);
        geoshow(latt,lont,'DisplayType','polygon','FaceColor',...
            [0.9,0.9,0.9],'EdgeColor','k')
        geoshow(scLat,scLong,'Color','black')
        grid on
    end

    
end

function [tab, solution] = generateTable(solutions,params)
    [~,ind] = sort(solutions(1,:));
    sortedSolns = solutions(:,ind);
    goodSolns = sortedSolns(:,sortedSolns(1,:)>100);
    
    varNames = {'Altitude [km]', 'Period [min]', 'Orbit Repeat #',...
        'Time Until Repeat [hrs]'};
    tableArray = zeros(length(goodSolns),4);
    
    inc = params.inc;
    R = params.R;
    mu = params.mu;
    
    tableArray(:,1) = goodSolns(1,:)';
    tableArray(:,2) = 2.*pi.*sqrt((goodSolns(1,:)+R).^3 ./ mu)./60;
    tableArray(:,3) = goodSolns(3,:);
    tableArray(:,4) = (tableArray(:,2).*tableArray(:,3))./60;
    
%     for i = 1:length(goodSolns)
%         tableArray(i,1) = goodSoln
%         T = 2*pi*sqrt((goodSolns(1,i)+R)^3 / mu)/60;
%         
%     end
    
    tab = array2table(tableArray,'VariableNames',varNames); 
    solution = [];
    solution.a = tableArray(2,1) + R;
    solution.alt = tableArray(2,1);
    solution.T = tableArray(2,2);
    solution.n = tableArray(2,3);
    solution.inc = inc;
    solution.Tn = tableArray(2,4);
    
end