clear variables;
%% NON DIMENSIONALIZATION 
global S h1 h2 l1 per dep; 
U = 10; H = 1e-6; L = 1e-3; mu = 4e-3; S = 1;

%% TEXTURES AND PAD
padpos = 1.0; per = 20; l1 = 0.2; dep = 1;
h1 = padpos; h2 = padpos + dep;

%% MESH DEFINITIONS

ax = -0.1; bx = 1.1;
nx = 2048; courant = 1.0;
xx = linspace(ax, bx, nx); dx = xx(2)-xx(1);
dt = courant*dx/(S/2);

%% SOLVER

SOLVER = 1; %% 1 GAUSS-SEIDEL, 0 JACOBI

%% MORE SIMULATION PARAMETERS

tol = 1e-8; ome_p = 1;
hgap = 500; 

p0 = 0.00;

T0 = 0.0-l1; Tfinal = 0.77; NT = floor(abs(Tfinal-T0)/dt);

T = linspace(T0,Tfinal,NT); dt = T(2) - T(1);

fluidleft = 20; %%capa inicial en magnitud de H

%% 

vlimx = [ax bx];
PLOTP = 1; fP = 1e2;

%% SETTING GEOMETRY

z = xx*0+padpos; z(xx<0) = hgap; z(xx>1.0) = hgap;

nxper = round(per/dx);
if nxper <= nx
    htex = repmat(ht(T0, xx(1:nxper)),[ 1 ceil(nx/nxper)]);
else
    htex = ht(T0, linspace(ax,ax+per,nxper));
end

hh = z + htex(1:nx); hhm = hh;

%% INITIAL CONDITIONS

pn = xx'*0+p0;

%%

titlesize = 14;

s = [0.5*(hh(1:end-1).^3+hh(2:end).^3) 777];
    
pk = pn; pkm = pn;
PP = zeros(NT,nx);

abort = false;
button1 = waitbar(0,'Please wait...','CreateCancelBtn','abort=true;');

t = T0;
it = 1;

pressre = zeros(NT,1); pressea = zeros(NT,1); 
%presseaode = zeros(NT,1); betaode = zeros(NT,1)+l1/2; 
beta = zeros(NT,1)+l1/2;
betavel = NaN+zeros(NT,1);

%for ode45
%Fbtt = @(t,bet) S/2*( 1 + (S*t*h2^3) ./ ( (h2^3-h1^3)*S*t+bet*h1^3 ) );

for it = 2:NT

    htex = circshift(htex,[1 2*courant]);
    hh = z + htex(1:nx);
    
    x1n = xx(find(htex,1))-dx/2;
    
    t = x1n/S;
    x1m = x1n-dx;
    
    kt = 2*dx^2/dt;
    
    pk = circshift(pn,[1 2*courant]); pkm = circshift(pn,[1 2*courant]);
    
    if t > 0
        
        %% HEUN
        betavelE = Fb(x1n , beta(it-1) + dt*Fb(x1m,beta(it-1))); % FORWARD EULER
        betavel(it) = mean([ Fb(x1m,beta(it-1)), betavelE ]); % HEUN
        beta(it) = beta(it-1) + dt*betavel(it);
        
        %% ODE45
%         [aux, betaaux] = ode45(Fbtt,[0 t], l1/2 );
%         betaode(it) = betaaux(end);
        %%
        
        IB = find(pn(xx<beta(it)+0.1)>1e-11,1,'last'); 
        
        pressea(it) = (beta(it)-x1n)*(2*Fb(x1n , beta(it))-S)*(h2-h1)/h2^3;
%         presseaode(it) = (betaode(it)-x1n)*(2*Fb(x1n , betaode(it))-S)*(h2-h1)/h2^3;
        
        pressre(it) = x1n*(h2-h1)/h1^3*(l1*h1^3/(x1n*h2^3+l1*h1^3));
        
    end
    
%     if mod(it,PLOTP) == 0, dibujar_re_ea_analitic; end
    
    waitbar(it/NT,button1,'working...');if abort;break; end
end

dibujar_re_ea_analitic;
delete(button1)