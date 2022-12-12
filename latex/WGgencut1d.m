% WGgen : generateur d'une grille reguliere de type 
% differences/volumes finis
% Auteurs JC Toussaint & L Bastard - BE EM

function g=WGgen
clc
clear all
close all

h=0.1;
g=uniform_grid(-5.0, 5.0, h);
g.lambda=1.;

% abscisses des extremites
% decrivant le guide
region = [-1  +1];
 
indice=3.5;
g=insert(g, region, indice);
end
