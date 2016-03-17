%% Programa que resuelve las ecuaciones de Navier-Stokes en un tubo de seccion 
%% variable con perfil de velocidad conocida en la entrada.
%%Octubre 2013, Alfredo Jaramillo. ICMC, USP.

clear variables;
num_fig = 501;
%% geometria y parametros del problema

L = 1;
H = 1; UL = 0; UH = 10;

mu = 1; rho = 1.0; % viscocidad, densidad
fx = 0.0; fy = 0.0; %componentes de la fuerza aplicada sobre el fluido
Uref = UL + UH;

Re = rho*Uref*H/mu;
    
UL = UL/Uref; UH = UH/Uref;

tfinal = 1000;
tolU = 1e-6; tolP=1e-9;

%% definicion de las mallas

NY = 30; NX = NY * L/H; %numero de volumenes en X e Y

xx = linspace(0, L, NX+1); yy = linspace(0, H, NY+1);
dx  = xx(2) - xx(1);

[YYp, XXp] = meshgrid(yy(1:end-1)+dx/2, xx(1:end-1)+dx/2); % malla presiones

%% parametros metodo
dt = 0.3*min([0.25*dx^2*Re 2/Uref^2/Re]) ; % dt condicionado para estabilidad

%% caso malla geometria del proyecto
II = reshape(1:NX*NY, NX, NY); % indices de todos los volumenes
IW = II(1,2:end-1)'; IE = II(end,2:end-1)'; 
IN = II(2:end-1,end); IS = II(2:end-1,1); %indices, oeste, este, norte, sur
INE = II(end,end); INW = II(1,end); ISE = II(end,1); ISW = II(1,1);

%%

II([IW;IE;IN;IS;INE;INW;ISE;ISW]) = -1; % marco los nodos no internos

II = II(:); II(II==-1) = []; %quito los nodos no internos

%% condiciones iniciales y de borde
Vn = 0*YYp; Un = 0*YYp; % velocidades verticales, horizontales (here is the no-split condition)
Pn = 0*XXp; % campo de presiones 

Utemp = Un; Vtemp = Vn; Ax = Utemp*0; Dx = Ax; Ay = 0*Vtemp; Dy = Ay;
Vnp = Vn; Unp = Un;
Pnp = Pn;

%% iteraciones, metodo explicito

crear_matriz_presiones3

niter = 0;
resU = Inf; resP = Inf;

salidas1 = []; salidas2=[];


while (niter*dt < tfinal && resU>tolU && resP >tolP)
%while 1
    resuelve_utemp_explicito3
    resuelve_presion
    resuelve_unp_3
    
    resU = (norm(Un-Unp,'inf')+norm(Vn-Vnp,'inf')) /...
           ( max(norm(Un,'inf'),UH) + max(norm(Vn,'inf'),UL) );
    %resP = norm(Pnp-Pn,'inf')/norm(Pn,'inf');
    
    if(mod(niter,60)==0)
        graficas
    end
    
    Un = Unp; Vn = Vnp; Pn = Pnp;
 
    niter=niter+1;
end

graficas

% grafico_salida

%     if(~mod(niter,258) || niter ==0 || niter==2584)
%         salidas1 = [salidas1;[niter Un(nxp1,:)]];
%         salidas2 = [salidas2;[niter Un(nxp2,:)]];
%     end 