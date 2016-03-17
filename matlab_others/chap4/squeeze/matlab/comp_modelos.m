se = load('sigma_ea_900.dat');
sr = load('sigma_re_900.dat');

hold on
plot(se(:,1),se(:,2),'--b')
plot(sr(:,1),sr(:,2),'--m')

Tp= linspace(0.25001,0.49999,400);
sig1 = 1-sqrt(p0*h(Tp).^3./dh(Tp));
plot(Tp,sig1,'-r');