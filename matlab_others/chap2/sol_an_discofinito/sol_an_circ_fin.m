R = 128; 
x0 = -0.5;
h0 = 0.5;
S = 1;
mu = 1.3e-2;

L = 1e-3; H = 1e-6;

pref = 6*mu*S*L/H^2;

one_atm = 101325;

gm0 = atan(x0/sqrt(2*R*h0*H/L));

xx = linspace(x0,-x0,400); dx = xx(2)-xx(1);
gg = atan(xx/sqrt(2*R*h0*H/L));

%%
pright = 0*38.5*one_atm/pref;
pleft = 0*one_atm/pref;

%% solucion analitica
f = @(gm) (pleft-pright)*h0^2/S/sqrt(2*R*H/L*h0) + 0.5*(gm-gm0) + 0.25*(sin(2*gm)-sin(2*gm0)) - 1./cos(gm).^2.*...
    ( 3/8*(gm-gm0) + .25*(sin(2*gm)-sin(2*gm0)) + 1/32*(sin(4*gm)-sin(4*gm0)) ) ;

gbar = fzero(f,-gm0/2);

C = pleft - S/h0^2*sqrt(2*R*H/L*h0) * ( 0.5*gm0 + 0.25*sin(2*gm0) - 1/cos(gbar)^2*( 3/8*gm0 + .25*sin(2*gm0) + 1/32*sin(4*gm0)) );

p = @(g) S/h0^2*sqrt(2*R*H/L*h0)*(0.5*g + 0.25*sin(2*g) - 1./cos(gbar)^2.*( 3/8*g + .25*sin(2*g) + 1/32*sin(4*g) )) + C;

%%

h = h0 + xx.^2/(2*R*H/L);
% 
% hold on
% 
[~,II] = min(abs(gg-gbar));
hflux = h(II); theta = xx; theta(1:II) = 1; theta(II+1:end) = hflux./h(II+1:end);
% 

pres = [p(gg(1:II)) 0*xx(II+1:end)+pright];
% 
gradp = [pres(2)/dx (pres(3:end)-pres(1:end-2))/dx/2 (pres(end)-pres(end-1))/dx];
pois = -h.^3/2.*gradp; coue = S/2*h.*theta;

% 
% figure(1)
% plot(xx, pres*pref/one_atm,'r--',xx,h);grid on

[ax,h1,h2]=plotyy(xx,pres,xx,h);

%xlabel('movement direction'); 
%ylim(ax(2),[0 max(h)])
%ylabel(ax(1),'hydrodynamic pressure');ylabel(ax(2),'pad profile');
%set(gcf,'color','w');


 figure(2)
plot(xx, pois, 'b--',xx,h,'k-', xx, coue,'r--', xx, pois+coue,'g--')
% grid on

% figure(3)
% hold on
% plot(xx,h-1,'b--','linewidth',2')
% plot(xx,R*L/H-sqrt(R^2-xx.^2)*L/H,'k-','linewidth',2')
% set(gcf,'color','w');set(gca,'FontSize',14)
% xlabel('movement direction [mm]'),ylabel('ring profile [\mu m]');
% grid on