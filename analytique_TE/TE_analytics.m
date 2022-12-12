function TE_analytics
clc
clear all
close all

um=1; % mum
lambda=1.0*um;
global a
a=0.3*um;

global k0
k0=2*pi/lambda;

global n1
global n2
global n3
n1=1.0; n2=3.5; n3=1.0;

neff=linspace(n1, n2, 1000);
plot(neff, f(neff)/pi, 'r-');
grid on

m=0;
lst_neff=[];
while true
    fun = @(neff) f(neff)-m*pi;
    try
        sol=fzero(fun, [n1, n2]);
        lst_neff=[lst_neff sol];
    catch me
        break
    end
    m=m+1;
end

lst_neff

figure
for i=1:length(lst_neff)
    neff=lst_neff(i);
    k1=k0*sqrt(neff^2-n1^2);
    k2=k0*sqrt(n2^2-neff^2);
    k3=k0*sqrt(neff^2-n3^2);
    phi=k2*a-atan(k3/k2);
    
    xtab=[];
    Etab=[];
    
    x=linspace(-2, -a, 10000);
    xtab=[xtab, x];
    E=cos(k2*a-phi)*exp(k3*(a+x));
    Etab=[Etab, E];
    
    x=linspace(-a, +a, 10000);
    xtab=[xtab, x];
    E=cos(k2*x+phi);
    Etab=[Etab, E];
    
    x=linspace(+a, +2, 10000);
    xtab=[xtab, x];
    E=cos(k2*a+phi)*exp(k1*(a-x));
    Etab=[Etab, E];
    
    subplot(length(lst_neff), 1, i) 
    plot(xtab, Etab);
    legend(num2str(neff))
    hold on;
    grid on
end

end

function [res]=f(neff)
global a
global k0
global n1
global n2
global n3

k1=k0*sqrt(neff.^2-n1^2);
k2=k0*sqrt(n2^2-neff.^2);
k3=k0*sqrt(neff.^2-n3^2);
res=2*k2*a-atan(k3./k2)-atan(k1./k2);  % m=0
end
