x = 0:1e-2:1;

r1 = 0.6*x' + 0.05;

p1 = (4*(1-x).*x)';

%% gauss seidel
nx = length(x); dx = x(2)-x(1);
DS = 1/dx^2*(diag(-2*ones(1,nx-2)) + diag(ones(nx-3,1),1));
In = 1/dx^2*diag(ones(nx-3,1),-1);
b = -8*ones(nx-2,1);

v = DS \ b;

solaux = p1*0.5;
err = 1e-5+1;
u = zeros(nx,1); 

while (err>1e-6)
    u(2:end-1) = -DS\In*solaux(2:end-1) + v;
    u( u > r1 ) = r1( u > r1 );
    err = norm(solaux-u,'inf')/norm(u,'inf');
    solaux = u;
end

plot(x,p1,x,r1)
hold on 
plot(x,u,'r')

axis image

