nt = 2626;
nx = 450;

x = linspace(0,1,nx);
t = linspace(0.245,2.245,nt);

p = zeros(nt,nx);
sigma = zeros(nt,1) + 0.5;

for i = 1:nt    
    p(i,:) = load(strcat('res',num2str(i,'%04d'),'.dat'));
    I = find(p(i,:)==0);
    if ~isempty(I)
        sigma(i) = x(I(end));
    else
        sigma(i) = 0.5;
    end
end

sigma(sigma < 0.5) = 0.5;



% for i=1:nt
%     plot(x,p(i,:))
%     drawnow
%     system('sleep 0.2');
% end

plot(t,sigma);
xlim([0.25 1.75])