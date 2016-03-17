%% solves for a membrane displacement with an obstacle

x = linspace(0,1,100);

r1 = 0*x' + 0.0;

p1 = (4*(1-x).*x)';

%% gauss seidel
nx = length(x); dx = x(2)-x(1);

DS = 1/dx^2*(diag(-2*ones(1,nx-2)) + diag(ones(nx-3,1),1));
F = 1/dx^2*diag(ones(nx-3,1),-1);

RHS = (x.*(1-x))'.*sin(5*pi*x)' +0.1;

v = DS \ RHS(2:end-1);

uwo = (DS + F ) \ RHS(2:end-1);

err = 1e-5+1;
u = zeros(nx,1); 

solaux = uwo;

while (err>1e-6)
    u(2:end-1) = -DS\F*solaux + v;
    u( u < r1 ) = r1( u < r1 );
    err = norm(solaux-u(2:end-1),'inf');
    solaux = u(2:end-1);
end

hold on 
plot(x,RHS)
plot(x,1e2*u,'--r','Linewidth',2)
plot(x,1e2*[0; uwo; 0],'-r','Linewidth',2)

set(gcf,'Color','w'); set(gca,'FontSize',14)
grid on

