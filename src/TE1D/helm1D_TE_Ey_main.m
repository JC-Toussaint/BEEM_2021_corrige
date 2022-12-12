% Stern Progress in Electromagnetics Research (1995), PIER 10, 123-186
% Finite Difference Analysis of Planar Optical Waveguides
%
% Quasi-TE (Ey) modes
%
%   ._____________.____._______.____.__________________.
%
% condtions de dirichlet

% helm1D_TE_Ey_main
clc
clear all
close all

h=0.001;
g=uniform_grid(-1.0, 1.0, h);
g.lambda=1;

a=0.3;
nG=3.5;

region = [-a +a];
g=insert(g, region, nG);

nmode=5;
g=TE_solve(g, nmode);

for n=1:nmode
    display_mode(g, n)
end


