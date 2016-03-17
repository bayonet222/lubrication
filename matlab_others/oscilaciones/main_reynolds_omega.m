clear variables;
%% NON DIMENSIONALIZATION 
global S h1 h2 l1 per dep; 
U = 10; H = 1e-6; L = 1e-3; mu = 4e-3; S = 1;

%% TEXTURES AND PAD
padpos = 1.0; per = 2.4; l1 = 0.2; dep = 1.0;
h1 = padpos; h2 = padpos + dep;

%% MESH DEFINITIONS

ax = 0.0; bx = 1.5;
nx = 2^8; courant = 1.0;
xx = linspace(ax, bx, nx); dx = xx(2)-xx(1);
dt = courant*dx/(S/2);

%% SOLVER

SOLVER = 1; %% 1 GAUSS-SEIDEL, 0 JACOBI

%% MORE SIMULATION PARAMETERS

tol = 1e-5; ome_p = 1;
hgap = 500; 

T0 = -l1; Tfinal = 1.0-l1;

T = T0:dt:Tfinal; NT = length(T);

fluidleft = 20; %%INITIAL FLUID HEIGHT

kt = 2*dx^2/dt; 

%% 

vlimx = [0 1];
PLOTP = 1; fP = 1e2;

%% SETTING GEOMETRY

z = xx*0+padpos; z(xx<0) = hgap; z(xx>1.0) = hgap;

htex = ht(T(1),xx);
hh = z + htex; hhm = hh;

%% INITIAL CONDITIONS

pn = xx'*0;

%%

titlesize = 14;
s = [0.5*(hh(1:end-1).^3+hh(2:end).^3) 777];
    
pk = pn; pkm = pn;
PP = zeros(NT,nx);

abort = false;
button1 = waitbar(0,'Please wait...','CreateCancelBtn','abort=true;');

beta = l1/2;
vbeta = []; betat=[]; vbetan = [];

t = T(1);
it = 1; dibujar_r
I = 2:nx-1; 

for it = 2:NT

    t = T(it);
    dt  = T(it)-T(it-1);
    
    kt = 2*dx^2/dt; 
    
    htex = ht(t, xx);
    hh = z + htex;
    
    s = [0.5*(hh(1:end-1).^3+hh(2:end).^3) 777];
    thetan(1) = (fluidleft + htex(1))/hh(1);
    
    pk = pn; pkm = pn;
    
    change = Inf;
    
    D =  sparse([1 I nx], [1 I nx], [1 -(s(I)+s(I-1)) 1], nx, nx);
    L =  sparse(I, I-1, s(I-1), nx, nx);
    U =  sparse(I, I+1, s(I), nx, nx);
    %b = [pk(1) S*dx*(hh(I)-hh(I-1)) - 2*S*(htex(I+1)-2*htex(I)+htex(I-1)) pk(end)]';
    b = [pk(1) S*dx*(hhm(I)-hhm(I-1)) + kt*(hh(I)-hhm(I)) pk(end)]';
    
    if SOLVER
        y = (D+L)\b; M = -(D+L)\U; %GAUSS-SEIDEL
    else
        y = D\b; M = -D\(L+U); %JACOBI
    end
    
    iter = 0;
    while (change > tol)
        iter = iter + 1;
        %% P
            % EXPLICITO
            
        pk = ome_p*(M*pk + y) + (1-ome_p)*pkm;
        
        pk(pk<0) = 0.0; % CAVITATION CHECK

        if mod(iter, 50) ==0
            change = norm(pk-pkm,'inf')/norm(pkm,'inf');
        end
             
        pkm = pk;
    end
    
    hhm = hh;
    pn = pkm;
    
    PP(it,:) = pn;
    
    if t > 0

        x1m = S*T(it-1); x1n = S*T(it);
        
        betatN = Fb(x1n , beta + dt*Fb(x1m,beta));
        betat = [betat; 0.5*( Fb(x1m, beta) + betatN)]; % HEUN
        %betat = Fb(x1,beta); % EULER
        
        beta = beta + dt*betat(end);

        vbeta = [vbeta; beta];
        
        IB = find(pn(xx<beta+0.1)>1e-11,1,'last'); vbetan = [vbetan; xx(IB)];
        
        pa = (beta-x1n)*(2*Fb(x1n , beta)-S)*(h2-h1)/h2^3;
        
%         betatode45 = Fbt(t,betaode45(tip));
%         paode45 = (betaode45(tip)-x1n)*(2*betatode45-S)*(h2-h1)/h2^3;
        
        %pb = 2*(S-betat(end))*(h2-h1)/h1^3*x1;
    end
    
    if mod(it,PLOTP)==0, dibujar_r; end
    
    waitbar(it/NT,button1,'working...');if abort;break; end
end

dibujar_r;
delete(button1)