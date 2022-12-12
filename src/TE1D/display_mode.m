function display_mode(g, nmode)
Nx=g.Nx;
X=g.x;

sol=g.V(:, nmode);

figure
plot(g.x, sol);
title(['neff=' num2str(g.neff(nmode))])
grid on
end

