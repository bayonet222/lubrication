%% VERSION SIN MODELO P-THETA, SOLO MULTIGRID PARA MODELO DE REYNOLDS


clear variables;
tic

tol = 1e-5; M = 7; 
NMg = 3; % numero de mallas

nx = 2.^(M:-1:M-NMg+1) + 1;
%% NON DIMENSIONALIZATION 
% H = 1e-6; L = 1; mu = 4e-3;
p0 = 0.025;

%% MESH DEFINITIONS

xx = linspace(0, 1.0, nx(1)); dx = 1./(nx-1);
dt = .3*dx(1);

%% upper surface

h0 = .375; h1 = .125; w = 4*pi;
H = @(t) (h1*cos(w*t) + h0);

%% MORE SIMULATION PARAMETERS
ome_p = 1;

mgItD = [1 1 1 1 30];

T0 = 0.0; Tfinal = 1.5  ; NT = floor((Tfinal-T0)/dt);
T = linspace(T0, Tfinal, NT); 

PLOTP = 2; fP = 1e0;
IRES = 5;

%% INITIAL CONDITIONS

pn = p0 + xx'*0;

%%    
pk = cell(NMg,1); b = pk;

for img = 1:NMg
    pk{img} = zeros(nx(img),1);
    b{img} = zeros(nx(img),1);
end

pk{1} = p0 + 0*pk{1};
pn = pk{1}; pkm = pn;

abort = false;

it = 1; 
I = (2:nx-1)'; 

h = H(0);

dibujar_squeeze

sigma = 0.5 + 0*T;

%A = sparse([1 2:nx(1)-1 nx(1) 2:nx(1)-1 2:nx(1)-1], ...
%                     [1 2:nx(1)-1 nx(1) 1:nx(1)-2 3:nx(1)], ...
%                     [1 -1/dx(1)^2*2*ones(1,nx(1)-2) 1 1/dx(1)^2*ones(1,nx(1)-2) 1/dx(1)^2*ones(1,nx(1)-2)]);full(A)

for it = 2:NT

    t = T(it);
    
    h = H(t); s = h^3;
        
    change = Inf;
    
    b{1} = [p0; 2/s*( -h1*w*sin(w*t) )*ones(nx(1)-2,1); p0];
    
    iter = 0;
    
    pk{1} = p0 + pk{1}*0;
    while (change > tol)

        iter = iter + 1;
        %% MULTIGRID DOWN
        
        pk{1} = gsC(b{1}, dx(1), 1, pk{1}, 1 );
        pk{1}(pk{1}<0)=0;
        for img = 2:NMg
            b{img} = interp_h_H( residueC(b{img-1}, pk{img-1}, dx(img-1)) );
            pk{img} = gsC(b{img}, dx(img), mgItD(img), pk{img}, 0);
            pk{img}(pk{img}<0)=0;
        end
        
        for img = NMg:-1:2
            pk{img-1} = pk{img-1} + interp_H_h( pk{img} );
        end

%         if (mod(iter,IRES)==0)
%             change = norm( residueC(b{1}, pk{1}, dx(1)), 'inf');
%         end
        change = norm(pk{1}-pkm,'inf');
        pkm = pk{1};
    end
    
    pn = pkm;
    
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

fprintf('\n');
dibujar_squeeze