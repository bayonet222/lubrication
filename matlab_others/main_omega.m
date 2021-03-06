clear variables;
%% NON DIMENSIONALIZATION 
global S h1 h2 l1 per dep; 
U = 10; H = 1e-6; L = 1e-3; mu = 4e-3; S = 1;

%% TEXTURES AND PAD
padpos = 1.0; per = 2.4; l1 = 0.2; dep = 1.0;
h1 = padpos; h2 = padpos + dep;

%% MESH AND DOMAIN DEFINITIONS

ax = 0; bx = 1.5;
nx = 2^9; courant = 0.9;
xx = linspace(ax, bx, nx); dx = xx(2)-xx(1);
dt = courant*dx/(S/2);

%% MORE SIMULATION PARAMETERS
tol = 1e-9; tolR = 1e-7;
ome_t = 1.0; ome_p = 1.0;
hgap = 500; 

T0 = -l1; Tfinal = 0.5; NT = floor((Tfinal-T0)/dt);
T = linspace(T0,Tfinal,NT);

%T = [T dt:dt:Tfinal]; NT = length(T);

fluidleft = 20; %%capa inicial en magnitud de H 

%% PLOT
vlimx = [ax bx];
PLOTP = 1; fP = 1e2;

%% SETTING GEOMETRY

z = xx*0+padpos; z(xx<0) = hgap; z(xx>1.0) = hgap;

htex = ht(T(1),xx);
hhn = z + htex;

%% INITIAL CONDITIONS

pn = xx*0; thetan = (fluidleft+htex)./hhn; thetan(thetan>1.0) = 1.0;

%%

pk = pn; thetak = thetan;% pkm = pn; thetakm = thetak;
cm = hhn.*thetan; ck = cm;

TH = zeros(NT,nx);
PP = zeros(NT,nx);

abort = false;
% button1 = waitbar(0,'Please wait...','CreateCancelBtn','abort=true;');

beta = l1/2;
vbeta = []; betat=[]; vbetan = [];

it = 1; x1 = 0;

%solving beta with ode45
%[~, betaode45] = ode45(@Fbt, [0 T(T>0)], l1/2);

tip = 1; 
II = 2:nx-1;

for it = 2:NT

    t = T(it);
    
    dt = T(it)-T(it-1); 
    kt = 2*dx^2/dt;
    
    %% VARIABLES EN EL TIEMPO N
    
    htex = ht(t, xx); %htex in time n
    hhn = z + htex;
    
    sn = [0.5*(hhn(1:end-1).^3+hhn(2:end).^3) 777];
        
    thetak(1) = 1.0;
    
    %%
    change = Inf; res = Inf;
    iter = 0;
    while ( (res > tolR) || isnan(res))
        %% P
        iter = iter+1;
      %  if (t>dt), pk = sol_e(xx,beta,t*S,Fb(S*t,beta),h1,h2); end
        
        IA = find((pk(2:end-1) > 0) | (thetak(2:end-1) >=1)) + 1;
        
        %PK = -kt*( ck(IA)-cm(IA) ) - S*dx*( ck(IA)-ck(IA-1) ) ... %IMPLICITO
        %             + sn(IA).*pk(IA+1) + sn(IA-1).*pk(IA-1) ;
        PK = -kt*( ck(IA)-cm(IA) ) - S*dx*( cm(IA)-cm(IA-1) ) ... %EXPLICITO
                 + sn(IA).*pk(IA+1) + sn(IA-1).*pk(IA-1) ;
        

        PK = PK./(sn(IA)+sn(IA-1));
        
        pk(IA) = ome_p*PK + (1-ome_p)*pk(IA);
        
        thetak( IA(pk(IA)>=0) ) = 1.0; pk(IA(pk(IA)<0)) = 0; % cavitation check
        
        %% THETA
        IC = find((pk(2:end-1)<=0) | (thetak(2:end-1)<1)) + 1;% IC([1 end])=[];
        
        % IMPLICITO
        %TK = kt*cm(IC) + S*dx*ck(IC-1) + sn(IC).*(pk(IC+1)-pk(IC))... 
        %     -sn(IC-1).*(pk(IC)-pk(IC-1));
        %TK = TK ./ ( (kt+S*dx)*hhn(IC) );
        
        % EXPLICITO
        TK = kt*cm(IC) - S*dx*(cm(IC)-cm(IC-1)) + sn(IC).*(pk(IC+1)-pk(IC))...
             - sn(IC-1).*(pk(IC)-pk(IC-1));
        TK = TK ./ (kt*hhn(IC));
        
        % RELAXATION
        thetak(IC) = ome_t*TK + (1-ome_t)*thetak(IC);
        
        %% ULTIMO NODO DERECHA
        thetak(end) = kt*cm(end) + S*dx*(cm(end-1)-cm(end)) + sn(end)*(0-pk(end))... % EXPLICITO
             -sn(end-1)*(0-pk(end-1)); % EXPLICITO
        thetak(end) = thetak(end) ./ (kt*hhn(end)); %EXPLICITO
        
        %thetak(end) = kt*cm(end) + S*dx*ck(end-1) + sn(end)*(0-pk(end))... %IMPLICITO
        %     -sn(end-1)*(0-pk(end-1));
        %thetak(end) = thetak(end) ./ ((dx*S+kt)*hhn(end)); %IMPLICITO
        
        %%
        
        pk(IC(thetak(IC)<1)) = 0; thetak(IC(thetak(IC)>=1)) = 1; % cavitation check
        

        if mod(iter,100) == 0
            LHS = kt*(ck(II)-cm(II)) + S*dx*(cm(II)-cm(II-1)); %EXPLICITO
            %LHS = kt*(ck(II)-cm(II)) + S*dx*(ck(II)-ck(II-1)); %IMPLICITO
            RHS = sn(II).*pk(II+1) - (sn(II)+sn(II-1)).*pk(II) + sn(II-1).*pk(II-1);  
            res = norm(LHS-RHS,'inf');
        end
        
        
        %change = (norm(pk-pkm,'inf') + norm(thetak-thetakm,'inf'))/...
        %         (norm(pkm,'inf') + norm(thetakm,'inf'));
             
        pkm = pk; thetakm = thetak; ck = hhn.*thetak;
        
    end  % while (change > tol)

    %% ACTUALIZACION DE VARIABLES AUXILIARES
    
    thetan = thetak; pn = pk;
    cm = ck;
    
    %%
    
    TH(it,:) = thetan; PP(it,:)=pk;
    
    if t > 0
        x1m = S*T(it-1); x1n = S*T(it);
        
        betatN = Fb(x1n , beta + dt*Fb(x1m,beta));
        betat = [betat; 0.5*( Fb(x1m, beta) + betatN)]; % HEUN
        
        beta = beta + dt*betat(end);

        vbeta = [vbeta; beta];
        
        IB = find(pn(xx<beta+0.1)>1e-11,1,'last'); vbetan = [vbetan; xx(IB)];
        
        pa = (beta-x1n)*(2*Fb(x1n , beta)-S)*(h2-h1)/h2^3;
        
%         betatode45 = Fbt(t,betaode45(tip));
%         paode45 = (betaode45(tip)-x1n)*(2*betatode45-S)*(h2-h1)/h2^3;
        
        %pb = 2*(S-betat(end))*(h2-h1)/h1^3*x1;
    end
    
    if mod(it,PLOTP)==0, dibujar; end

%     waitbar(it/NT, button1,'working...');
    if abort,
         break; 
    end
end

dibujar;
% delete(button1)