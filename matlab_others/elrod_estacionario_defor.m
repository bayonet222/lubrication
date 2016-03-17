p0 = 0.0;

nx =  512;
x = linspace(0,1,nx)'; dx = x(2)-x(1);

%C = -sum(1./h.^2) / sum(1./h.^3);

A1 = .25; A2 = 0.16;
h = 1.25 + ( (1-x)*A1 + A2*x ).*cos(2*2*pi*x) + (A1-A2)*x;
%h = 1.25 + .25.*cos(2*pi*x);

%%
% h = 2*x;
% h(1:round(nx/3)) = h(round(nx/3));
% h(round(2*nx/3):end) = h(round(2*nx/3));
% %%

theta = ones(nx,1);

s = [.5*(h(1:end-1).^3 + h(2:end).^3); h(end)^3];

c = h;

p = zeros(nx,1); p(1) = p0; p(end)=p0;

tol = 1e-7;
err = tol + 1;

auxp = p; auxt = theta;
while (err>tol)
    
    for i=2:nx-1  
        if (p(i) > 0 || theta(i) == 1.0)
            p(i) = 1/(s(i)+s(i-1))*(-dx*(c(i)-c(i-1)) + s(i)*auxp(i+1) + s(i-1)*auxp(i-1) );
            if (p(i)>=0)
                theta(i) = 1.0;
            end
        end
        
        if(p(i)<=0 || theta(i)<1.0)
            theta(i) = 1/(dx*h(i))*(s(i)*(auxp(i+1)-auxp(i)) - s(i-1)*(auxp(i)-auxp(i-1)) + dx*c(i-1));    
            if (theta(i)<1.0)
                p(i) = 0.0;
            else
                theta(i) = 1.0;
            end
        end
    end
    
    err = norm(auxp-p,'inf')+norm(auxt-theta,'inf');
    auxp = p;
    auxt = theta;
    c = theta.*h;
end

J = -h.^3.*[diff(p); p(1)-p(end)]/dx + h.*theta;
% p = ( DpS + In ) \ b;

figure(2)
plot(x,p/max(p),x,h)
