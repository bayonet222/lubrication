global S h1 h2 l1 per dep; 
U = 10; H = 1e-6; L = 1e-3; mu = 4e-3; S = 1;
%% TEXTURES
per = 2.4; l1 = 1/10; dep = 1;

%% MESH DEFINITIONS
nx = 512; courant = 0.5;
xx = linspace(0,2,nx); dx = xx(2)-xx(1);
dt = dx*courant;

%% NUMERIC PARAMETERS AND SIMULATION
tol = 1e-5; ome_t = 1.0; ome_p = 1.0;
hgap = 500; S = 1;

NT = 2000;
padpos = 1;
theta0 = (padpos+0.5)/hgap;
PLOTP = 1;

%% INITIAL CONDITIONS
pn = xx*0; thetan = xx*0+theta0;


%%
kt = 2*dx^2/dt; 
t = 0; 

z = xx*0+padpos; z(1:floor(nx/3)) = hgap; z(floor(nx/3*2):end) = hgap;
thetan(z==padpos) = 1.0;

htex = ht(0,xx);
hh = z+htex;
cn = hh.*thetan;

s = [0.5*(hh(1:end-1).^3+hh(2:end).^3) 777];
    
pk = pn; thetak = thetan; pkm = pn; thetakm = thetak;

for it = 1:NT

    t = t + dt;
    htex = ht(t,xx);
    hh = z + htex;
    s = [0.5*(hh(1:end-1).^3+hh(2:end).^3) 777];
    
    pk = pn; thetak = thetan; pkm = pn; thetakm = thetak;
    ck = hh.*thetan;
    
    change = Inf;
    
    while (change > tol)
        %% pressure
        
        %IA = find((pk > 0) | (thetak==1));
        IA = find((pk > 0) | (thetak>=1));
        
        PK = -kt*( ck(IA)-cn(IA) ) - S*dx*( ck(IA)-ck(IA-1) ) ...
                 + s(IA).*pk(IA+1) + s(IA-1).*pk(IA-1) ;
        PK = PK./(s(IA)+s(IA-1));
        
        pk(IA) = ome_p*PK + (1-ome_p)*pk(IA);
        
        thetak( IA(pk(IA)>=0) ) = 1.0; pk(IA(pk(IA)<0)) = 0; % cavitation check
        %thetak(pk>=0) = 1.0; pk(pk<0) = 0; % cavitation check
        
        %% theta
        IC = find((pk<=0) | (thetak<1)); IC([1 end])=[];
        
        TK = kt*cn(IC) + S*dx*ck(IC-1) + s(IC).*(pk(IC+1)-pk(IC))...
             -s(IC-1).*(pk(IC)-pk(IC-1));
        TK = TK ./ ((kt+S*dx).*hh(IC));
        
        thetak(IC) = ome_t*TK + (1-ome_t)*thetak(IC);
        
        thetak(end) = kt*cn(end) + S*dx*ck(end-1) + s(end)*(0-pk(end))...
             -s(end-1)*(0-pk(end-1));
        thetak(end) = thetak(end) ./ ((kt+S*dx).*hh(end));

        pk(IC(thetak(IC)<1)) = 0; thetak(IC(thetak(IC)>=1)) = 1; % cavitation check
        
        change = (norm(pk-pkm,'inf') + norm(thetak-thetakm,'inf'))/...
                 (norm(pkm,'inf') + norm(thetakm,'inf'));
        pkm = pk; thetakm = thetak; ck = hh.*thetak;
    end
    thetan = thetak; pn = pk;
    cn = hh.*thetan;
    
    if mod(it,PLOTP)==0, dibujar; end

end