% encuentra las velocidades del paso n+1

dPU = diff([Pnp;Pnp(1,:)]);
Unp(IPU)  = Utemp(IPU) - dt/dx * dPU(IPU-1);  %nodos internos, pagina 29, ecuacion (2.19)
Unp([IW;INW;ISW]) = Utemp([IW;INW;ISW]) - dt/dx * dPU(end,:)';

dPV = diff(Pnp,1,2);
Vnp(IPV)  = Vtemp(IPV) - dt/dx * dPV(IPV);