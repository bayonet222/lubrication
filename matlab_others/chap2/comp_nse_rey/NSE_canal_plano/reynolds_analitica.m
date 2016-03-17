x = 0:0.01:1;

S = 1; 
L = 1e-4;% [0.1*mm] %largo del dominio periódico
H = 1e-6;% [micrómetro]

h0 = 10; %h min
h1 = 1; %profundidad de la perturbación

h = @(x) h0 - h1/2*(1-cos(2*pi*x));

%% condición p(0)=p(L)=0
%% debemos determinar c1 de la ecuación
%% p(L) = \int_0^1 { (c1+S*h) / h^3 } dx = 0


% c1 = -9.465:0.0001:-9.46;
% root = c1*0;
% for i=1:length(c1)
%     root(i) = integral(@(x) h0*(c1(i)+S*h(x))./h(x).^3,0,1);
% end
% plot(c1,root,'*')


c1 = -9.4606;
f = @(x) h0*(c1+S*h(x))./h(x).^3;

p = 0;
for i=1:length(x),
    p(i)=integral(f,0,x(i));
end;

plot(x,p*1e3,x,-(h(x)-h0)); grid on
