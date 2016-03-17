%clear variables;
tic

tol = 1e-5; nx = 900;
%% NON DIMENSIONALIZATION 
% H = 1e-6; L = 1; mu = 4e-3;
p0 = 0.025; S = 0.0;

%% MESH DEFINITIONS

xx = linspace(0, 1.0, nx); dx = xx(2)-xx(1);
dt = .3*dx;

%% upper surface

h0 = .375; h1 = .125; w = 4*pi;
H = @(t) (h1*cos(w*t) + h0)*ones(nx,1);
dH = @(t) (-w*h1*sin(w*t))*ones(nx,1);

h = @(t) (h1*cos(w*t) + h0);
dh = @(t) (-w*h1*sin(w*t));
%H = @(t) 0.01*t*ones(nx,1)+0.3;

%% MORE SIMULATION PARAMETERS
ome_p = 1;

T0 = 0.0; Tfinal = 1.0; NT = floor((Tfinal-T0)/dt);
T = linspace(T0, Tfinal, NT); 

kt = 2*dx^2/dt; 

PLOTP = 1; fP = 1e0;

%% INITIAL CONDITIONS

pn = xx'*0+p0;

%%    

pk = pn; pkm = pn;
PP = zeros(NT,nx);

abort = false;

it = 1; 
I = (2:nx-1)'; 

hm = H(0); hh = H(0);

dibujar_squeeze

sigma = 0.5 + 0*T;

for it = 2:NT

    t = T(it);
    
    hh = H(t);
    
    s = hh(I).^3;
    
    pk = pn*0+p0;
    
    change = Inf;
    
    TL =  sparse([1; I; nx; I], [1; I; nx; I-1], [1; -2*s; 1; s], nx, nx);
    U =  sparse(I, I+1, s, nx, nx);
    b = [p0; kt*(hh(I)-hm(I)); p0];
    
    pkm = (TL+U)\b;
    
    if ( pkm(end/2)<=0 ) %% por simetría si p(end/2)>=0 no es necesario GS
        y = TL \ b; M = - TL \ U;
        pkm(pkm<0) = 0;
        while (change > tol)
            paux = M*pkm + y;

            pk = ome_p*paux + (1-ome_p)*pk;

            pk(pk<0) = 0; % CAVITATION CHECK

            change = norm(pk-pkm,'inf')/norm(pk,'inf');

            pkm = pk;
        end
    end
    
    hm = hh;
    pn = pkm;
    PP(it,:) = pn;
    
%     if t>0
%         x1=S*t;
%         betat = [betat; ...
%                  0.5*( Fb(x1, beta) + ...
%                  Fb(x1+dt*S, beta + dt*Fb(x1,beta)) )]; % HEUN
%         %betat = Fb(x1,beta); % EULER
%         
%         beta = beta + dt*betat(end); vbeta = [vbeta; beta];
%         IB = find(pn(xx<beta+0.1)>1e-11,1,'last'); vbetan = [vbetan; xx(IB)];
%         
%         pa = (beta-x1)*(2*betat(end)-S)*(h2-h1)/h2^3;
%         pb = 2*(S-betat(end))*(h2-h1)/h1^3*x1;
%     end
    
    if mod(it,PLOTP)==0, dibujar_squeeze; end
    
    k = mod(it/NT,0.1);
    if k >= 0.1-1/NT;
        fprintf('%1.0f%% ',it/NT*100);
    end
    
    isigma = find(pn==0,1,'last');
    if ~isempty(isigma)
        sigma(it) = xx(isigma);
    else
        sigma(it) = 0.5;
    end
end
toc

% fprintf('\n');
dibujar_squeeze