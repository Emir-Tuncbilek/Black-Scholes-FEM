function [kQ,cQ,bQ,fQ]=timeEvals(xq, t)
sigma = 0.2;
r = 0.05;
Smax = 50;                         % from main.m  ==> change it for real production code !!!!

kQ = 0.5 * sigma.^2 * xq.^2;       % diffusion coefficient
cQ = (sigma^2 - r) * xq;           % convection coefficient
bQ = r;                            % reaction term

%====================================%
%======= Manufactured Solution ======%
%====================================%

u = (1-t) * xq * (Smax - xq);
ut = -xq * (Smax - xq);
us = (1-t) * (Smax - 2*xq);
uss = -2 * (1-t);

%====================================%
% Manufactured Solution placed on fQ %
%====================================%

fQ = 0; % ut -( (.5 * sigma^2 * xq^2 * uss) + (r * xq * us) - (r * u));

end