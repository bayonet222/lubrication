x0 = 0;
h0 = 1.0;
S = 1;
mu = 1.3e-2;

L = 1e-3; H = 1e-6;
R = 80*H/L; 

pref = 6*mu*S*L/H^2;

one_atm = 101325;

K = sqrt(2*R*h0);


xx = linspace(-0.5,0.5  ,400); dx = xx(2)-xx(1);
gg = atan(xx/K);

A = [-S*K*(3/8*gg(1)+sin(2*gg(1))/4 + sin(4*gg(1))/32) 1;
    -S*K*(3/8*gg(end)+sin(2*gg(end))/4 + sin(4*gg(end))/32) 1];
b = [-S*K*(gg(1)/2+sin(2*gg(1))/4);-S*K*(gg(end)/2+sin(2*gg(end))/4)]; 

aux = A\b;

C = aux(2); gbar = acos(sqrt(1/aux(1)));

%%

p = @(g) S*K*(0.5*g + 0.25*sin(2*g) - 1./cos(gbar)^2.*( 3/8*g + .25*sin(2*g) + 1/32*sin(4*g) ));
%%

h = h0 + xx.^2/(2*R);
% 
% hold on
% 
[~,II] = min(abs(gg-gbar));
hflux = h(II); theta = xx; theta(1:II) = 1; theta(II+1:end) = hflux./h(II+1:end);
% 

pres = p(gg);
% 
% gradp = [pres(2)/dx (pres(3:end)-pres(1:end-2))/dx/2 (pres(end)-pres(end-1))/dx];
% pois = -h.^3/2.*gradp; coue = S/2*h.*theta;

% 
% figure(1)
% plot(xx, pres*pref/one_atm,'r--',xx,h);grid on

[ax,h1,h2]=plotyy(xx,pres,xx,h);

xlabel('movement direction'); 
ylim(ax(2),[0 max(h)])
ylabel(ax(1),'hydrodynamic pressure');ylabel(ax(2),'pad profile');
set(gcf,'color','w');
grid on

% figure(2)
% plot(xx, pois, 'b--',xx,h,'k-', xx, coue,'r--', xx, pois+coue,'g--')
% grid on

% figure(3)
% hold on
% plot(xx,h-1,'b--','linewidth',2')
% plot(xx,R*L/H-sqrt(R^2-xx.^2)*L/H,'k-','linewidth',2')
% set(gcf,'color','w');set(gca,'FontSize',14)
% xlabel('movement direction [mm]'),ylabel('ring profile [\mu m]');
% grid on