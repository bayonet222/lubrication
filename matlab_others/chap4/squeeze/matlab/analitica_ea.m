se = load('sigma_ea_4096.dat');
sr = load('sigma_re_900.dat');

hold on
%plot(se(:,1),se(:,2)-0.5,'--b')
%plot(sr(:,1),sr(:,2),'--m')

% tref = fzero(@(t)3*h1*sin(w*t)^2+(h1*cos(w*t)+h0)*cos(w*t),0.31);
tref = 0.314825861279561; TFIN = 0.7325;

p0 = 0.025;
h0 = .375; h1 = .125; w = 4*pi;

h = @(t) h1*cos(w*t) + h0;
dh = @(t) -w*h1*sin(w*t);

SIG = @(t) 1-sqrt(p0*h(t).^3./dh(t)) ;

%% CREA MALLA PARA INTERPOLACION de H(SIGMA)

t = linspace(0.250076,tref,2000);
HH = h(t); SS = SIG(t);

hIN = @(s) interp1(SS,HH,s);

%%

FSIG = @(t,s) 0.5*( (1-s).^2.*dh(t) - p0*h(t).^3 ) ./ ( (s-1).*(hIN(s)-h(t)) ) ;

Tp = linspace(0.250079,0.499365,1000);
plot(Tp,SIG(Tp),'-r'); ylim([0.5 1]); xlim([0.2 0.76])

tref = 0.31483;
TINT = linspace(tref,TFIN,500);
[aux,sigNUM] = ode23(FSIG,TINT,SIG(tref));

plot(aux,sigNUM,'-b');