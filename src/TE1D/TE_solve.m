function g=TE_solve(g, n)
%
% g=TE_solve(g, nmode)
% calcul des n valeurs propres neff
% et leur mode associe
%
% entree : structure g et nb de modes a calculer
% sortie : structure g contenant les donnees membres
% vecteur des valeurs propres D
% vecteur des neff
% tableaux des modes propres ranges par colonnes

g = dirichlet(g);
A = build_A(g);

N=g.N;
ld=g.ld;

Id=speye(N ,N);

P=Id;
P(ld, :)=[];

%RESOLUTION DE L'EQUATION AUX VALEURS PROPRES
indice_max=max(g.indice(:));
A=A-indice_max^2*Id;  % decalage de valeur propre
Ap=+P*A*P';

[Vp, D] = eigs(Ap, n, 'sm');
% vecteurs propres associes
g.V=P'*Vp;

% valeurs propres recorrigees
g.D=diag(D)+indice_max^2;
g.neff=sqrt(g.D);
end

function g = dirichlet(g)
ld = [1 g.Nx];
g.ld=unique(ld);
end


