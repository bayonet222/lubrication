%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% nodos internos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pagina 30 Prosperetti & Tryggvason
Ax(II) = .25/dx*( (Un(II)+Un(II+1)).^2 - (Un(II-1)+Un(II)).^2 ...
            + (Un(II+NX)+Un(II)) .* (Vn(II-1)+Vn(II)) ...
            - (Un(II-NX)+Un(II)) .* (Vn(II-1-NX)+Vn(II-NX)) );
        
Ay(II) = .25/dx*( (Un(II+NX+1)+Un(II+1)) .* (Vn(II)+Vn(II+1)) ...
         - (Un(II+NX)+Un(II)) .* (Vn(II)+Vn(II-1))...
         + (Vn(II+NX)+Vn(II)).^2 - (Vn(II)+Vn(II-NX)).^2 );

% Pagina 31

Dx(II) = 1/dx^2*(Un(II+1)+Un(II-1)+Un(II+NX)+Un(II-NX)-4*Un(II));
Dy(II) = 1/dx^2*(Vn(II+1)+Vn(II-1)+Vn(II+NX)+Vn(II-NX)-4*Vn(II));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% velocidades horizontales borde
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sur Un's
Ax(US) = .25/dx*( (Un(US)+Un(US+1)).^2 - (Un(US-1)+Un(US)).^2 ...
            + (Un(US+NX)+Un(US)) .* (Vn(US-1)+Vn(US)) ...
            );
            %- (Un(US-NX)+Un(US)) .* (Vn(US-1-NX)+Vn(US-NX)) ); % no
            %hay flujo al sur

Dx(US) = 1/dx^2*(Un(US+1)+Un(US-1)+a*(8/3)+(4/3)*Un(US+NX)-6*Un(US));

% norte Un's
Ax(UN) = .25/dx*( (Un(UN)+Un(UN+1)).^2 - (Un(UN-1)+Un(UN)).^2 ...
            + 0 ...% no hay flujo al norte
            - (Un(UN-NX)+Un(UN)) .* (Vn(UN-1-NX)+Vn(UN-NX)) ); 
Dx(UN) = 1/dx^2*(Un(UN+1)+Un(UN-1)+(4/3)*Un(UN-NX)-6*Un(UN));

% izquierda de este Un's IW
Ax(UE) = .25/dx*( (Un(UE)+UnE(UE/NX)).^2 - (Un(UE-1)+Un(UE)).^2 ...
            + (Un(UE+NX)+Un(UE)) .* (Vn(UE-1)+Vn(UE)) ...
            - (Un(UE-NX)+Un(UE)) .* (Vn(UE-1-NX)+Vn(UE-NX)) );
Dx(UE) = 1/dx^2*(UnE(UE/NX)+Un(UE-1)+Un(UE+NX)+Un(UE-NX)-4*Un(UE));
% 
% % Un NE
Ax(UNE) = .25/dx*( (Un(UNE)+UnE(UNE/NX)).^2 - (Un(UNE-1)+Un(UNE)).^2 ...
            + 0 ... % no hay flujos al norte
            - (Un(UNE-NX)+Un(UNE)) .* (Vn(UNE-1-NX)+Vn(UNE-NX)) );
Dx(UNE) = 1/dx^2*(UnE(UNE/NX)+Un(UNE-1)+(4/3)*Un(UNE-NX)-6*Un(UNE));
% 
% % Un SE
Ax(USE) = .25/dx*( (Un(USE)+UnE(USE/NX)).^2 - (Un(USE-1)+Un(USE)).^2 ...
            + (Un(USE+NX)+Un(USE)) .* (Vn(USE-1)+Vn(USE)) ...
            - 0); % no hay flujos al sur
Dx(USE) = 1/dx^2*(UnE(USE/NX)+Un(USE-1)+(8/3)*a+(4/3)*Un(USE+NX)-6*Un(USE));

%
% % IEw
Ax(IEw) = .25/dx*( 0 - (Un(IEw-1)+Un(IEw)).^2 ... %no hay flujo al este (pared)
            + (Un(IEw+NX)+Un(IEw)) .* (Vn(IEw-1)+Vn(IEw)) ...
            - (Un(IEw-NX)+Un(IEw)) .* (Vn(IEw-1-NX)+Vn(IEw-NX)) );
Dx(IEw) = 1/dx^2*(Un(IEw+NX)+Un(IEw-NX)+(4/3)*Un(IEw-1)-6*Un(IEw));
   
% % INEw
Ax(INEw) = .25/dx*( 0 - (Un(INEw-1)+Un(INEw)).^2 ... %no hay flujo al este (pared)
            + 0 ... %no hay flujo al norte (pared)
            - (Un(INEw-NX)+Un(INEw)) .* (Vn(INEw-1-NX)+Vn(INEw-NX)) );
        
Dx(INEw) = 1/dx^2*((4/3)*Un(INEw-NX)+(4/3)*Un(INEw-1)-8*Un(INEw));
% este Un's ( derivada normal de la velocidad al lado este igual a 0 )
% (pagina 36, Zero gradient condition)
UnEtemp = 4/3*Un([USE;UE;UNE]) - 1/3*Un([USE;UE;UNE]-1); % aproximacion de segundo orden de la ecuación dUn/dx = 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% velocidades verticales borde
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sur Vn's
Ay(VS) = .25/dx*( (Un(VS+NX+1)+Un(VS+1)) .* (Vn(VS)+Vn(VS+1)) ...
         - (Un(VS+NX)+Un(VS)) .* (Vn(VS)+Vn(VS-1))...
         + (Vn(VS+NX)+Vn(VS)).^2 - (Vn(VS)+0).^2 ); %VN(VS-NX)=0, velocidad vertical en las paredes
Dy(VS) = 1/dx^2*(Vn(VS+1)+Vn(VS-1)+Vn(VS+NX)+0-4*Vn(VS));

% Vn SUROESTE: velocidad vertical nula en el lado OESTE y SUR
Ay(VSW) = .25/dx*( (Un(VSW+NX+1)+Un(VSW+1)) .* (Vn(VSW)+Vn(VSW+1)) ...
         - 0 ...     
         + (Vn(VSW+NX)+Vn(VSW)).^2 - (Vn(VSW)+0).^2 );
     
Dy(VSW) = 1/dx^2*((4/3)*Vn(VSW+1)+Vn(VSW+NX)+0-6*Vn(VSW));

% Vn's OESTE: velocidad vertical nula en el lado OESTE
Ay(VW) = .25/dx*( (Un(VW+NX+1)+Un(VW+1)) .* (Vn(VW)+Vn(VW+1)) ...
         - 0 ... 
         + (Vn(VW+NX)+Vn(VW)).^2 - (Vn(VW)+Vn(VW-NX)).^2 );
     
Dy(VW) = 1/dx^2*((4/3)*Vn(VW+1)+Vn(VW+NX)+Vn(VW-NX)-6*Vn(VW));

% Pagina 29
Utemp = Un + dt*(-Ax + mu*Dx + fx);
Vtemp = Vn + dt*(-Ay + mu*Dy + fy);

% Vn ESTE y SURESTE
% aproximacion de segundo orden de la ecuación dVn/dx = 0
% Vtemp([VE; VSE]) = 3/2*Vn([VE;VSE]-1)-1/2*Vn([VE;VSE]-2); 

% % aproximacion de segundo orden de la ecuación V(ESTE) = 0
Vtemp([VE; VSE]) = 1/(15*dx^2)*(6*Vn([VE;VSE]-2)-10*Vn([VE;VSE]-1))*(2*dx)^2 ...
                   - 1/(15*dx) *(21*Vn([VE;VSE]-2)-25*Vn([VE;VSE]-1))*(2*dx)+Vn([VE;VSE]-2);

% Vtemp([VE; VSE]) = 0;
