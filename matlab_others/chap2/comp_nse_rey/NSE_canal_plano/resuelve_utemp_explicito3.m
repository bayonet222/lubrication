%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NODOS INTERNOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ADVECCION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pagina 30 Prosperetti & Tryggvason
Ax(II) = .25/dx*( (Un(II) + Un(II+1)).^2 - (Un(II-1) + Un(II)).^2 ...
            + (Un(II+NX) + Un(II)) .* (Vn(II-1) + Vn(II)) ...
            - (Un(II-NX) + Un(II)) .* (Vn(II-1-NX) + Vn(II-NX)) );
        
Ax(IS) = .25/dx*( (Un(IS) + Un(IS+1)).^2 - (Un(IS-1) + Un(IS)).^2 ...
            + (Un(IS+NX) + Un(IS)) .* (Vn(IS-1) + Vn(IS)) ...
            ); % - (Un(IS-NX) + Un(IS)) .* (Vn(IS-1-NX) + Vn(IS-NX)) ); % no hay flujo al sur        
            
Ax(IN) = .25/dx*( (Un(IN) + Un(IN+1)).^2 - (Un(IN-1) + Un(IN)).^2 ...
            + 0 ...% no hay flujo al norte
            - (Un(IN-NX) + Un(IN)) .* (Vn(IN-1-NX) + Vn(IN-NX)) ); 

Ax(IW) = .25/dx*( (Un(IW) + Un(IW+1)).^2 - (Un(IE) + Un(IW)).^2 ... %PERIOCIDAD -> IW - 1 = IE
            + (Un(IW+NX) + Un(IW)) .* (Vn(IE) + Vn(IW)) ...
            - (Un(IW-NX) + Un(IW)) .* (Vn(IE-NX) + Vn(IW-NX)) );

% izquierda de este Un's IW
Ax(IE) = .25/dx*( (Un(IE) + Un(IW)).^2 - (Un(IE-1) + Un(IE)).^2 ... %PERIOCIDAD -> IE + 1 = IW
            + (Un(IE+NX) + Un(IE)) .* (Vn(IE-1) + Vn(IE)) ...
            - (Un(IE-NX) + Un(IE)) .* (Vn(IE-1-NX) + Vn(IE-NX)) );

% % Un NE
Ax(INE) = .25/dx*( (Un(INE) + Un(INW)).^2 - (Un(INE-1) + Un(INE)).^2 ... %PERIOCIDAD -> INE + 1 = INW
            + 0 ... % no hay flujos al norte
            - (Un(INE-NX) + Un(INE)) .* (Vn(INE-1-NX) + Vn(INE-NX)) );  

% % Un SE
Ax(ISE) = .25/dx*( (Un(ISE) + Un(ISW)).^2 - (Un(ISE-1) + Un(ISE)).^2 ...%PERIOCIDAD -> ISE + 1 = ISW
            + (Un(ISE+NX) + Un(ISE)) .* (Vn(ISE-1) + Vn(ISE)) ...
            - 0); % no hay flujos al sur
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ay(II) = .25/dx*( (Un(II+NX+1) + Un(II+1)) .* (Vn(II) + Vn(II+1)) ...
         - (Un(II+NX) + Un(II)) .* (Vn(II) + Vn(II-1))...
         + (Vn(II+NX) + Vn(II)).^2 - (Vn(II) + Vn(II-NX)).^2 );
     
% sur Vn's
Ay(IS) = .25/dx*( (Un(IS+NX+1) + Un(IS+1)) .* (Vn(IS) + Vn(IS+1)) ...
         - (Un(IS+NX) + Un(IS)) .* (Vn(IS) + Vn(IS-1))...
         + (Vn(IS+NX) + Vn(IS)).^2 - (Vn(IS) + 0).^2 ); %IN(IS-NX)=0, velocidad vertical en las paredes

% Vn's OESTE
Ay(IW) = .25/dx*( (Un(IW+1+NX) + Un(IW+1)) .* (Vn(IW) + Vn(IW+1)) ...
         - (Un(IW+NX) + Un(IW)) .* (Vn(IW) + Vn(IE))...
         + (Vn(IW+NX) + Vn(IW)).^2 - (Vn(IW) + Vn(IW-NX)).^2 );
     
% Vn SUROESTE
Ay(ISW) = .25/dx*( (Un(ISW+NX+1) + Un(ISW+1)) .* (Vn(ISW) + Vn(ISW+1)) ...
         - (Un(ISW+NX) + Un(ISW)) .* (Vn(ISW) + Vn(ISE))... 
         + (Vn(ISW+NX) + Vn(ISW)).^2 - (Vn(ISW) + 0).^2 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DIFUSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pagina 31
Dx(II) = 1/dx^2*( Un(II+1) + Un(II-1) + Un(II+NX) + Un(II-NX) - 4*Un(II) );

Dx(IS) = 1/dx^2*(Un(IS+1) + Un(IS-1) + 4/3*Un(IS+NX) + 8/3*UL - 6*Un(IS));
Dx(IN) = 1/dx^2*(Un(IN+1) + Un(IN-1) + 4/3*Un(IN-NX) + 8/3*UH - 6*Un(IN));

Dx(IW) = 1/dx^2*( Un(IW+1) + Un(IE) + Un(IW+NX) + Un(IW-NX) - 4*Un(IW) );
Dx(ISW) = 1/dx^2*( Un(ISW+1) + Un(ISE) + 4/3*Un(ISW+NX) + 8/3*UL - 6*Un(ISW) );
Dx(INW) = 1/dx^2*( Un(INW+1) + Un(INE) + 4/3*Un(INW-NX) + 8/3*UH - 6*Un(INW) );

Dx(IE) = 1/dx^2*( Un(IW) + Un(IE-1) + Un(IE+NX) + Un(IE-NX) - 4*Un(IE) );
Dx(INE) = 1/dx^2*(Un(INW) + Un(INE-1) + 4/3*Un(INE-NX) + 8/3*UH - 6*Un(INE));
Dx(ISE) = 1/dx^2*(Un(ISE) + Un(ISE-1) + 4/3*Un(ISE+NX) + 8/3*UL - 6*Un(ISE));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dy(II) = 1/dx^2*( Vn(II+1) + Vn(II-1) + Vn(II+NX) + Vn(II-NX) - 4*Vn(II) );

Dy(IS) = 1/dx^2*( Vn(IS+1) + Vn(IS-1) + Vn(IS+NX) + 0 - 4*Vn(IS) );

Dy(IW) = 1/dx^2*( Vn(IW+1) + Vn(IE) + Vn(IW+NX) + Vn(IW-NX) - 4*Vn(IW));
Dy(ISW) = 1/dx^2*( Vn(ISW+1) + Vn(ISE) + Vn(ISW+NX) + 0 - 4*Vn(ISW) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pagina 29
Utemp = Un + dt*( -Ax + 1/Re*Dx + fx );
Vtemp = Vn + dt*( -Ay + 1/Re*Dy + fy );