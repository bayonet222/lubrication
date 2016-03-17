
x0 = -0.5;
h0 = 1;
S = 1;
mu = 1.3e-2;

L = 1e-3; H = 1e-6;

pref = 6*mu*S*L/H^2;

one_atm = 101325;

xx = linspace(x0,-x0,100); dx = xx(2)-xx(1);


%%

%Rs = [2 4 6 8 16:16:256];
h0s = [0.5:0.5:10];


prights = zeros(length(h0s),length(Rs));
Ropt = zeros(1,length(h0s));

for ih0 = 1:length(h0s)
    R = 2;
    diffR = Inf;
    prighta = 0;
    
    while (diffR > 0)
        R = R + 1;
        
        h0 = h0s(ih0);

        gm0 = atan(x0/sqrt(2*R*h0*H/L));

        pright = 0;
        dpright = 5*one_atm/pref;
        tolE = 1e-3*one_atm/pref;
        
        while (dpright > tolE)

        %% solucion analitica
            f = @(gm) -pright*h0^2/S/sqrt(2*R*H/L*h0) + 0.5*(gm-gm0) + 0.25*(sin(2*gm)-sin(2*gm0)) - 1./cos(gm).^2.*...
                ( 3/8*(gm-gm0) + .25*(sin(2*gm)-sin(2*gm0)) + 1/32*(sin(4*gm)-sin(4*gm0)) );

            gbar = fzero(f,-gm0/2);            
            
            if  (0 < gbar) && (gbar < 0.5)
                pright = pright + dpright;
            else
                pright = pright - dpright;
                dpright = dpright/2;
            end
        end

        diffR = pright - prighta;
        prighta = pright;
    end
    
    Ropt(ih0) = R - 1;
    disp(ih0/length(h0s))
end

plot(h0s,Ropt,'x--','Linewidth',2,'MarkerSize',10); grid on

xlabel('h_0','FontSize',14);ylabel('R_{opt}','FontSize',14)

set(gca,'FontSize',14)
set(gcf,'color','w'); grid on