function [sigma] = fixSigma(a)
%fixSigma Summary of this function goes here
%   Detailed explanation goes here
    Pk = 2*pi*sqrt(a^3/mu);
    Pomega = Pk*(1-1.5*J2*(R/a)^2 * (3 - 4*sind(i)^2));
    p = a*(1-e^2);
    omegaDot = -1.5*(n*R^2 * J2)/(p^2) * cosd(i);
    lambdaDot = omegaDot - w;
    k = linspace(0,100);
    sigma = (k.*Pomega.*lambdaDot)./(2*pi);
end