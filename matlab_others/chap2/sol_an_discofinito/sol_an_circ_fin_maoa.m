
x0 = -0.5;
h0 = 1;
S = 1;
mu = 1.3e-2;

L = 1e-3; H = 1e-6;

pref = 6*mu*S*L/H^2;

one_atm = 101325;

xx = linspace(x0,-x0,100); dx = xx(2)-xx(1);


%%

Rs = [2 4 6 8 16:16:256]; h0s = [0.5 1.0 1.5 2.0 5 7.5 10.0];
Rs = 50; h0s = 0.2:0.2:10;


prights = zeros(length(h0s),length(Rs));

for ih0 = 1:length(h0s)
    for iR = 1:length(Rs)

        h0 = h0s(ih0);
        R = Rs(iR);

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

        prights(ih0,iR) = pright*pref/one_atm;
    end
    ih0/length(h0s)
end

plot(Rs,log10(prights)','x-','Linewidth',2);

xlim(Rs([1 end]));
xlabel('R','FontSize',14);ylabel('log(p_{max})','FontSize',14)

set(gca,'FontSize',14)
set(gcf,'color','w'); grid on

M=[];for i = 1:length(h0s);M{i}=strcat('h_0=',num2str(h0s(i)));end
% [RR, HH0] = meshgrid(Rs, h0s);
% surf(RR,HH0, log10(prights))
% surf(RR,HH0, prights)

% C = - S/h0^2*sqrt(2*R*H/L*h0) * ( 0.5*gm0 + 0.25*sin(2*gm0) - 1/cos(gbar)^2*( 3/8*gm0 + .25*sin(2*gm0) + 1/32*sin(4*gm0)) );

% p = @(g) S/h0^2*sqrt(2*R*H/L*h0)*(0.5*g + 0.25*sin(2*g) - 1./cos(gbar)^2.*( 3/8*g + .25*sin(2*g) + 1/32*sin(4*g) )) + C;

%%

% h = h0 + xx.^2/(2*R*H/L);
% 
% hold on
% 
% [~,II] = min(abs(gg-gbar));
% hflux = h(II); theta = xx; theta(1:II) = 1; theta(II+1:end) = hflux./h(II+1:end);

% pres = [p(gg(1:II)) 0*xx(II+1:end)+pright];

% gradp = [pres(2)/dx (pres(3:end)-pres(1:end-2))/dx/2 (pres(end)-pres(end-1))/dx];
% pois = -h.^3/2.*gradp; coue = S/2*h.*theta;





% figure(1)
% plot(xx, pres*pref/1000,'r--',xx,h);grid on
% 
% figure(2)
% plot(xx, pois, 'b--',xx,h,'k--', xx, coue,'r--', xx, pois+coue,'g--')
% grid on

