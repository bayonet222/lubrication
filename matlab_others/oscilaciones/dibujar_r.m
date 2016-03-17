clf
hold on
plot(xx, fP*pn,'r--');
plot(xx, -htex, 'k-','linewidth',2);
plot(xx, z, 'k-','linewidth',2)
ylim([-dep padpos+2]);
xlim(vlimx)

x1n = S*T(it);

if (T(it)>0)
    plot([0 x1n beta xx(end)],[0 fP*pa 0 0],'-b', 'linewidth',2);
    stit = sprintf('T %1.2e, PRe %1.2e, PE %1.2e',T(it), fP*max(pn), fP*pa);
    title(stit,'FontSize',titlesize);
end
drawnow