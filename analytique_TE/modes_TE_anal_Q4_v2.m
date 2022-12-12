% analytique TE guide symetrique

function lst_neff=modes_TE_anal_Q4_v2
clc
clear all
close all

um=1;
lambda = 1*um;
global a;
a = 0.3*um;
global n;
n = 3.5;

global k0;
k0=2*pi/lambda;

neff=linspace(1, n, 10000);
fneff = f(neff);
plot(neff, fneff/pi, 'r');
grid on

m=0;
lst_neff=[]
for m=0:10
    fun = @(neff) f(neff)-m*pi;
    try
    sol=fzero(fun, [1, n]);
    lst_neff=[lst_neff sol];
    catch me
        break
    end
end
end

function [res]=f(neff)
global a
global k0
global n
res=2*k0*a*sqrt(n^2-neff.^2)-2*atan(sqrt((neff.^2-1)./(n^2-neff.^2)));
end
