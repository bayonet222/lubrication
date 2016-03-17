clf
hold on
plot(xx, fP*pn,'r*-');
plot(xx, -htex, 'k-','linewidth',2);
plot(xx, z, 'k-','linewidth',2)
plot(xx, thetan.*hhn-htex,'b--','linewidth',2);

ylim([-dep padpos+1]);
xlim(vlimx)

x1n = S*T(it);
if (T(it)>dt)
    plot([0  x1n beta 1],[0 fP*pa 0 0],'r-', 'linewidth',2);%plot(x1+.5,1e2*pb,'*g');
%     plot([0 x1n betaode45(tip)],[0 fP*paode45 0],'b-', 'linewidth',2);%plot(x1+.5,1e2*pb,'*g');
end
drawnow