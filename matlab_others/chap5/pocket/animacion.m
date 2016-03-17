proffiles = dir(strcat(pwd,'/profilex*'));

xr = zeros(length(proffiles),1)+0.5;
dt = 3.91E-4;

%for i = 120:250
for i = 1:2:length(proffiles)
    hold on
    aux = importdata(strcat(pwd,'/',proffiles(i).name));
    plot(aux.data(:,1)/0.1,1e3*aux.data(:,2))
    plot(aux.data(:,1)/0.1,aux.data(:,4),aux.data(:,1),aux.data(:,5))
    plot(aux.data(:,1)/0.1,aux.data(:,3).*aux.data(:,4),'--r');
    I = find(aux.data(:,2)==0.0);
    
%     if ~isempty(I)
%         xr(i) = (aux.data(max(I),1)-0.5);
%     end
    
    ylim([-2 3]);
    
    drawnow
%     system('sleep 0.05');
%     print('-depsc',strcat(num2str(i),'.eps'))
    clf();
end
hold off
t = dt*[1:length(proffiles)];
% plot(t,xr)
