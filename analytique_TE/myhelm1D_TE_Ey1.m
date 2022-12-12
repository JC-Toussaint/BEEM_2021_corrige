% Stern Progress in Electromagnetics Research (1995), PIER 10, 123-186
% Finite Difference Analysis of Planar Optical Waveguides
%
% TE (Ey) modes
%
%   ._____________.____._______.____.__________________.
%
% condtions de dirichlet 

function g=myhelm1D_TE_Ey1
clc
clear all
close all

um=1;

h=0.01*um;
g=uniform_grid(-5.0, 5.0, h);
g.lambda=1.0*um;

W =0.3*um;
nG=3.5*um;

nodes = [-W W];
g=insert(g, nodes, nG);

nmode=4;
g=solve(g, nmode);

for n=1:nmode
    display(g, n)
end
end

function g=uniform_grid(xmin, xmax, hx)
g.x=xmin+hx/2:hx:xmax-hx/2;

g.p     = g.x(:);
g.Nx=length(g.x);
g.reg =zeros(g.Nx, 1);
g.indice =ones(g.Nx, 1);
g.N=g.Nx;
g.nrg=0;
g.dx=hx;

figure(1)
plot(g.p, g.indice, 'ro');
hold on
end

function g=insert(g, nodes, indice)
assert(length(nodes)==2, 'node length > 2')
figure(1)

xmin=nodes(1);
xmax=nodes(2);
xmin=min(xmin, xmax);
xmax=max(xmin, xmax);
lst_in=find((g.p>xmin) & (g.p<xmax));

g.nrg=g.nrg+1;
g.reg(lst_in)=g.nrg;
g.indice(lst_in)=indice;

plot(g.p, g.indice, 'b-')
plot(g.p(lst_in), ones(length(lst_in), 1), 'bo-')
grid on
hold on
title('indices')
drawnow()
end

function g = dirichlet(g)
ld = [index(g, 1) index(g, g.Nx)];
g.ld=unique(ld);
end

function g=solve(g, nmode)
g = dirichlet(g);
A = build_A(g);

N=g.N;
ld=g.ld;

Id=speye(N ,N);

P=Id;
P(ld, :)=[];

%RESOLUTION DE L'EQUATION AUX VALEURS PROPRES
indice_max=max(g.indice(:));
A=A-indice_max^2*Id;
Ap=+P*A*P';

[Vp, D] = eigs(Ap, nmode, 'sm');
% vecteur propre associee valeur propre numero n
%solp = Vp(:,nmode); 
% solp = Vp(:,1);
% sol=P'*solp;
% g.sol=sol;
g.V=P'*Vp;
g.D=diag(D)+indice_max^2;
g.neff=sqrt(g.D);
end

function A=build_A(g)
Nx=g.Nx;
N =g.N;
dx=g.dx;
ld=g.ld;
k0=2*pi/g.lambda;

A=sparse(N, N);

for ix=2:Nx-1
    n=index(g, ix);   
    if any(ld==n) 
        continue;
    end
    
    A(n, n)=-2/dx^2;
    A(n, n)=A(n, n)+g.indice(ix)^2*k0^2;
    
    p=index(g, ix+1);
    A(n, p)=+1/dx^2;  
    
    p=index(g, ix-1);
    A(n, p)=+1/dx^2;
          
end
A = A/k0^2;

end

function n=index(g, i)
n=i;
end

function display(g, nmode)
Nx=g.Nx;
X=g.x;

sol=g.V(:, nmode);

figure
plot(g.x, sol);
title(['neff=' num2str(g.neff(nmode))])
grid on
end

