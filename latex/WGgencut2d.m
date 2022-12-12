% WGgen : generateur d'une grille reguliere de type 
% differences/volumes finis
% Auteurs JC Toussaint & L Bastard - BE EM

% dependance : fonction inpoly
% Auteur : Darren Engwirda
%
% Bibliographie :
% Stern Progress in Electromagnetics Research (1995), 
% PIER 10, 123-186
% Finite Difference Analysis of Planar Optical Waveguides
%

function g=WGgen
clc
clear all
close all

h=0.1;
g=uniform_grid(-5.0, 5.0, -3.0, 3.0, h, h);
g.lambda=1.;

% coordonnees des sommets du polygone
% decrivant une portion du guide
region = [
       -2  -1;
       +2  -1;
       +1  +1;
       -1  +1
       ];
 
indice=3.5;
g=insert(g, region, indice);
end

function g=uniform_grid(xmin, xmax, ymin, ymax, hx, hy)
% generation d'un reseau orthorhombique de noeuds
% chacun place au centre d'une maille rectangulaire
% de cote g.dx=hx et g.dy=hy
% g.Nx noeuds selon Ox
% g.Ny noeuds selon Oy
% g.indice(ix, iy) indice de refraction au noeud
% de coordonnees entieres (ix, iy) initialise a 1

lx=xmin+hx/2:hx:xmax-hx/2;
ly=ymin+hy/2:hy:ymax-hy/2;

[g.x,g.y] = meshgrid(lx, ly);
g.x=transpose(g.x);
g.y=transpose(g.y);

g.p     = [g.x(:),g.y(:)];
[g.Nx, g.Ny]=size(g.x);

g.reg =zeros(g.Nx, g.Ny);
g.indice =ones(g.Nx, g.Ny);
size(g.indice)

g.N=g.Nx*g.Ny;
g.nrg=0;
g.dx=hx;
g.dy=hy;
figure(1)
plot(g.p(:, 1), g.p(:, 2), 'ro');
hold on
end

function g=insert(g, node, indice)
% insertion d'un polygone decrivant une portion 
% de guide en 2D.
% node   : tableau 2D des coordonnees des sommets du polygone
% indice : valeur de l'indice de refraction du materiau
% associe au polygone.

figure(1)
tic, [in,bnd] = inpoly(g.p,node); t1 = toc;
if find(bnd==1)
    warning('nodes on boundaries detected')
end;
g.nrg=g.nrg+1;
lst_in=find(in==1);
g.reg(lst_in)=g.nrg;
g.indice(lst_in)=indice;
plot(g.p(in,1),g.p(in,2),'bo',g.p(bnd,1),g.p(bnd,2),'go')
hold on
title('Inside (blue), outside (red) & boundary (green)')

% Plot geometry
nn    = size(node,1);
cnect = [(1:nn-1)',(2:nn)'; nn,1];
patch('faces',cnect,'vertices',node,'facecolor','none',...
   'edgecolor','k'), axis equal 
drawnow()
end
