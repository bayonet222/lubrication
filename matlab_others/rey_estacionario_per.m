p0 = 0.001;

nx =  1024;
x = linspace(0,1,nx)'; dx = x(2)-x(1);

%C = -sum(1./h.^2) / sum(1./h.^3);

%% sinusoides decrecientes
A1 = .25; A2 = 0.25;
h = 1.25 + ( (1-x)*A1 + A2*x ).*cos(2*2*pi*x) + (A1-A2)*x;

% %% semi-circulos
% h = x*0 + 1.25;
% per = 1;
% l = 0.3;
% for i=1:length(h)
%     aux = mod(x(i),per)/per;
%     if (aux > l && aux < (1-l))
%        h(i) = 1.25 + sqrt(- (aux-0.5)^2 + (0.5-l)^2);
%     end
% end

% sinusoides uniformes
% h = 1.25 + .25.*cos(4*2*pi*x);
%

%% triangulos
% h = x*0 + 1.25;
% per = 1/2;
% l = 0.3;
% for i=1:length(h)
%     aux = mod(x(i),per)/per;
%     if (aux > l && aux < (1-l))
%        h(i) = 1.25 + abs(aux-0.5) - 0.5 + l;
%     end
% end
% %%

s = .5*[(h(1:end-1).^3 + h(2:end).^3); h(end-1)^3+h(1)^3];

DpS = sparse(diag(s(2:end)+s(1:end-1)) - diag(s(2:end-1),1));
DpS(1,end) = -s(end-1);

In  = sparse(-diag(s(2:end-1),-1));
In(end,1) = -s(end-1);

b = -dx*diff([h(1:end-1);h(1)]);

p = zeros(nx-1,1);

tol = 1e-12;
err = tol + 1;

b(1)=p0; DpS(1,:) = [1 zeros(1,nx-2)];


v = DpS \ b;

while (err>tol)
    aux = p;
    p = v - DpS\(In*p);
    p(p<0) = 0;
    err = norm(aux-p,'inf')/norm(p,'inf');
end

%%
p = [p;p(1)];

J = -h(1:end-1).^3.*diff(p)/dx + h(1:end-1);
% p = ( DpS + In ) \ b;

figure(2)
plot(x,p/max(p),x,h,x(1:end-1),J)
