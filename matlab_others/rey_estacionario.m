p0 = 0.0;

nx =  8;
x = linspace(0,1,nx)'; dx = x(2)-x(1);

%C = -sum(1./h.^2) / sum(1./h.^3);

%% sinusoides decrecientes
A1 = .25; A2 = 0.16;
h = 1.25 + ( (1-x)*A1 + A2*x ).*cos(2*2*pi*x) + (A1-A2)*x;

%% semi-circulos
% h = x*0 + 1.25;
% per = 1;
% l = 0.2;
% for i=1:length(h)
%     aux = mod(x(i),per)/per;
%     if (aux > l && aux < (1-l))
%        h(i) = 1.25 + sqrt(- (aux-0.5)^2 + (0.5-l)^2);
%     end
% end

%% sinusoides uniformes
% h = 1.25 + .25.*cos(8*2*pi*x);
%%

% %% triangulos
% h = x*0 + 1.25;
% per = 1/2;
% l = 0.3;
% for i=1:length(h)
%     aux = mod(x(i),per)/per;
%     if (aux > l && aux < (1-l))
%        h(i) = 1.25 + abs(aux-0.5) - 0.5 + l;
%     end
% end
% 
% h = 2*x;
% h(1:round(nx/3)) = h(round(nx/3));
% h(round(2*nx/3):end) = h(round(2*nx/3));

%%

% %% semi-triangulos
% h = x*0 + 1.25;
% per = 1;
% l = 0.3;
% for i=1:length(h)
%     aux = mod(x(i),per)/per;
%     if (aux > l && aux < (1-l))
%        h(i) = 1.25 + 2*(aux-l);
%     end
% end
%

s = [.5*(h(1:end-1).^3 + h(2:end).^3); h(end)^3];

DpS = sparse(diag(s(2:end-1)+s(1:end-2))-diag(s(2:end-2),1));
In  = sparse(-diag(s(2:end-2),-1));

b = - dx *(h(2:end-1)-h(1:end-2)); 
b(1) = b(1) + s(1)*p0; b(end) = b(end) + s(end)*p0;

p = (DpS+In)\b;

tol = 1e-12;
err = tol + 1;

v = DpS \ b;

while (err>tol)
    aux = p;
    p = v - DpS\(In*p);
    p(p<0) = 0;
    err = norm(aux-p,'inf')/norm(p,'inf');
end


p = [p0; p ; p0];

J = -h(1:end-1).^3.*diff(p)/dx + h(1:end-1);
% p = ( DpS + In ) \ b;

figure(1)
plot(x,p*90,x,h,x(1:end-1),J)
