x0 = 0;
h0 = 1.0;
S = 1;
mu = 1.3e-2;

L = 1e-3; H = 1e-6;
R = 80*H/L; 

pref = 6*mu*S*L/H^2;

one_atm = 101325;

K = sqrt(2*R*h0);

hold on
err = [];
set(gcf,'Color','w'); 

expo = 6:9;
err = 0*expo; tauA = 0*expo; tauf = 0*expo; normAinv = 0*expo;

for k = 1:length(expo)
    N = 2^expo(k);
    xx = linspace(-0.5,0.5,N+1); dx = xx(2)-xx(1);

    h = h0 + xx.^2/(2*R);
    dh = S*xx/R;

    %% NUMERICAL SOLUTION

    a = 0.5*(h(2:end).^3+h(1:end-1).^3);
    A = N^2*(-diag(a(1:end-1)+a(2:end)) + diag(a(2:end-1),1) + diag(a(2:end-1),-1));
    f = S*diff(h(1:end-1)')/dx;

    presnum = A \ f;


    %% ANALYTICAL SOLUTION
    gg = atan(xx/K);

    B = [-S*K*(3/8*gg(1)+sin(2*gg(1))/4 + sin(4*gg(1))/32) 1;
        -S*K*(3/8*gg(end)+sin(2*gg(end))/4 + sin(4*gg(end))/32) 1];
    b = [-S*K*(gg(1)/2+sin(2*gg(1))/4);-S*K*(gg(end)/2+sin(2*gg(end))/4)]; 

    aux = B \ b;

    C = aux(2); gbar = acos(sqrt(1/aux(1)));

    p = @(g) S*K*(0.5*g + 0.25*sin(2*g) - 1./cos(gbar)^2.*( 3/8*g + .25*sin(2*g) + 1/32*sin(4*g) ));
    pres = p(gg);

    %% COMPARISON
    plot(xx+0.5,1e2*[0;presnum;0],'-r');
    
    err(k) = sqrt(dx)*norm(pres(2:end-1)'-presnum,2);
    normAinv(k) = 1/norm(A,2);
    tauf(k) = sqrt(dx)*norm( dh(2:end-1)' - diff(h(1:end-1)')/dx , 2);
    tauA(k) = sqrt(dx)*norm( dh(2:end-1)' - A*pres(2:end-1)', 2);
    
end
plot(xx+0.5,1e2*pres,'-b');
plot(xx+0.5,h,'k-');
ylabel('pressure profiles','Fontsize',14); grid on

% figure(2)
% plot(log2(1./(2.^expo)),log2(tauA),'-')