%% resuelve el problema de Poisson de la pagina 29, ecuacion (2.20)

Ub = diff([zeros(NX,1) Vtemp],1,2) + diff([Utemp; Utemp(end,:)]) ; % V = 0 en la pared inferior
Ub(INW) = 0; %P(INW)=0

Pnp(IP) = dt/dx * Ap(IP,IP) \  Ub(IP) ;