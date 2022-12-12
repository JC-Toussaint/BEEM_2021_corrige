function A=build_A(g)
Nx=g.Nx;
N =g.N;
dx=g.dx;
ld=g.ld;
k0=2*pi/g.lambda;
sq=(k0*dx)^2;

A=sparse(N, N);

for ix=1:Nx
% Numerotation des lignes de A = celle des noeuds en 1D
    n=ix;
    if any(ld==n)  % methode generale applicable au 2D
        continue;
    end;

    A(n, n)=-2/sq+g.indice(ix)^2;

    p=ix+1;
    A(n, p)=+1/sq;

    p=ix-1;
    A(n, p)=+1/sq;

end

end



