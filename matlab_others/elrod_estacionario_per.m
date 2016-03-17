p0 = 0.03;

nx =  512;
x = linspace(0,1,nx)'; dx = x(2)-x(1);

%C = -sum(1./h.^2) / sum(1./h.^3);

A1 = .25; A2 = 0.25;
h = 1.25 + ( (1-x)*A1 + A2*x ).*cos(2*2*pi*x) + (A1-A2)*x;
%h = 1.25 + .25.*cos(2*pi*x);

theta = ones(nx,1);

s = .5*[(h(1:end-1).^3 + h(2:end).^3); h(1)^3+h(2)^3];

c = h;

p = zeros(nx,1)+p0;% p(1) = p0; p(end) = p(1);

im = round(0.5/dx);

p(im) = p0;

tol = 1e-8;
err = tol + 1;

auxp = p; auxt = theta;


while (err>tol)
    
    p(1) = 1/(s(1)+s(nx-1))*(-dx*(c(1)-c(nx-1)) + s(1)*auxp(1+1) + s(nx-1)*auxp(nx-1) );
    
    if p(1) < 0
        p(1) = 0;
    end
    
  
    for i = [2:im-1 im+1:nx-2]
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

    p(nx-1) = 1/(s(nx-1)+s(nx-2))*(-dx*(c(nx-1)-c(nx-2)) + s(nx-1)*auxp(1) + s(nx-2)*auxp(nx-2) );
    
    if p(nx-1) <=0
        p(nx-1) = 0;
    end
    
    p(nx) = p(1);
    
    err = norm(auxp-p,'inf')+norm(auxt-theta,'inf');
    auxp = p;
    auxt = theta;
    c = theta.*h;
end

J = -h.^3.*[diff(p); p(2)-p(1)]/dx + h.*theta;
% p = ( DpS + In ) \ b;

figure(2)
plot(x,p/max(p),x,h,x,J)
