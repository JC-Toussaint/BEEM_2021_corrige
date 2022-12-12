% analytique TE guide symetrique
clc
clear all
close all

um=1;
lambda = 1*um;
a = 0.3*um;
n = 3.5;

k0=2*pi/lambda;
neff=linspace(1, n, 10000);
fneff = 2*k0*a*sqrt(n^2-neff.^2)-2*atan(sqrt((neff.^2-1)./(n^2-neff.^2)));
plot(neff, fneff/pi, 'r');
grid on

% neff_ana  neff_num
% 3.425     3.425
% 3.193     3.194
% 2.774     2.776
% 2.090     2.095