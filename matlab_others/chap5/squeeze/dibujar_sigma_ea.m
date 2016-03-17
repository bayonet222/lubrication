
figure(2)
hold on
set(gcf,'Color','w');

sigma = load('ea_450/sigma.dat');

% sigma(sigma==0.5)=NaN;

plot(sigma(:,1),sigma(:,3),'-k','LineWidth',2)

analitica_ea

xlim([0.22 0.75])

%xlabel('time'); ylabel('\Sigma(t)');