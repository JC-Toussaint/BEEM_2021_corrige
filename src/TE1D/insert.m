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

