hold on
% plot(xx, fP*pn,'r-*');
plot(xx, -htex(1:nx), 'k-','linewidth',2);
plot(xx, z, 'k-','linewidth',2)
ylim([-dep padpos+2]);
xlim(vlimx)

x1n = S*t;

if (t>0)
    plot([0 x1n x1n+l1 xx(end)],[0 fP*pressre(it) 0 0],'-b', 'linewidth',2);
    plot([0 x1n beta(it) xx(end)],[0 fP*pressea(it) 0 0],'-r', 'linewidth',2);
%     plot([0 x1n betaode(it) xx(end)],[0 fP*presseaode(it) 0 0],'-b', 'linewidth',2);
    %stit = sprintf('T %1.2e, PRe %1.2e, PE %1.2e',t, fP*max(pn), fP*pauxteo);
    %title(stit,'FontSize',titlesize);
    set(gcf,'Color','w');set(gca,'FontSize',14);
    xlabel('x','FontSize',14); ylabel('z','FontSize',14);
end
drawnow
grid on

hold off