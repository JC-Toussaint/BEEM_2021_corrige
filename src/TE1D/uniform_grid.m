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

