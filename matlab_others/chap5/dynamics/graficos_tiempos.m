times = [17 42 77];

for in = 1:1
    hold on;set(gcf,'Color','w'); set(gca,'FontSize',14)
    
    pres = load(strcat('ea0',num2str(times(in)),'.dat'));
    pres(:,1) = pres(:,1)*10-0.5;
    
    plot(pres(:,1),1e3*pres(:,2),'--r','LineWidth',2)
    plot(pres(:,1),pres(:,4),'-k','LineWidth',2)
    plot(pres(:,1),pres(:,5),'-k','LineWidth',2)
end

for in = 2:length(times)
    pres = load(strcat('ea0',num2str(times(in)),'.dat'));
    pres(:,1) = pres(:,1)*10-0.5;
    
    plot(pres(:,1),1e3*pres(:,2),'--r','LineWidth',2)
    plot(pres(:,1),pres(:,5),'-k','LineWidth',2)
end

grid on
xlim([-0.1 1.1]);ylim([-1 3])
xlabel('x');ylabel('z')

%% COMBINE WITH main_reynolds.m