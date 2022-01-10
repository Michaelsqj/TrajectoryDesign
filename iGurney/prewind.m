function [g,k,s] = prewind(ks,gs,GMAX,SMAX,TS);

g1n = ceil(gs/(SMAX*TS));
g1 = (g1n:-1:0)/g1n*gs;
g = g1;
kn = ks+sum(g1(2:end))*4258*TS;
g2n = ceil(sqrt(kn/SMAX/TS/4258/TS));
g3n = 0;
if (g2n*SMAX*TS>GMAX)
  g2n = floor(GMAX/SMAX/TS);
  kn2 = ks+sum(g1)*4258*TS-g2n*g2n*4258*TS*SMAX*TS;
  g3n = ceil(kn2/4258/GMAX/TS);
end
g2 = [(1:g2n) g2n*ones(1,g3n) ((g2n-1):-1:0)]*SMAX*TS;
g2 = g2*kn/(sum(g2)*4258*TS);
sum(g2)*4258*TS;

g = [g1 -g2];
k = ks+cumsum(g(2:end))*4258*TS;
s = 0;
