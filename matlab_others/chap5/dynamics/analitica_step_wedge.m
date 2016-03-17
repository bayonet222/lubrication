%% ANALITICAL SOLUTION OF REYNOLDS EQUATION FOR
%% THE FINITE STEP WEDGE

mu = 0.013; U = 10;
H = 1e-6; L = 1e-3;
S = 1;

pref = 6*mu*U*L/H^2;

% NX = 1000;
% x = linspace(0,1,NX);

% l = 4/(sqrt(27)+9); h1 = (sqrt(3)+2)/2;

l = 0.1:0.01:0.3; h1 = 1.5:0.01:2.5;

[l,h1] = meshgrid(l,h1);

pmax = @(l,h1) S*(h1-1).*(1-l).*l./(1+l.*(h1.^3-1));
% gradI = pmax/(1-l); gradD = pmax/-l;

% II0 = x < 1-l; II1 = x >= 1-l;

% p = x;
% p(II0) = x(II0)*S*gradI; p(II1) = (x(II1)-1+l)*S*gradD + (1-l)*S*gradI;

W = @(l,h1) 0.5 * pmax(l,h1);
F = @(l,h1) 3*S*(h1-1).^2.*l.*(1-l)./(1+l.*(h1.^3-1)) + S*((1-l)./h1+l);

Cf = F(l,h1)./W(l,h1);

figure(1)
pcolor(1./l-1,h1,pmax(l,h1));xlabel('l');ylabel('h1');
shading interp

figure(2)
pcolor(1./l-1,h1,Cf);xlabel('l');ylabel('h1');
shading interp
