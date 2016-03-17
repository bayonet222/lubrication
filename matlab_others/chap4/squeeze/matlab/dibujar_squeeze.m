
figure(2)
clf
hold on
plot(xx, fP*pn,'r-','linewidth',2);
% plot(xx, hh.*ones(nx(1),1), 'k-','linewidth',2);
%ylim([0 h1+h0*1.05]);
ylim([-0.5 1]);

% theta = 1+0*xx; theta(pn==0) = 0.0;
% plot(xx, theta'.*hh,'b--','linewidth',2);

drawnow