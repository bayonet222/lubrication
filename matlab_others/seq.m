DP = 1;

z = xx*0+padpos; z(1:floor(nx/4)) = hgap; z(floor(nx/4*3):end) = hgap;

abort = false;
h = waitbar(0,'Please wait...','CreateCancelBtn','abort=true;');

for i =1:DP:NT,
    hti = ht((i-1)*dt,xx,per,l1,dep,S);
    hhi = z + htex;
    clf
    hold on
    plot(xx, 1e1*PP(i,:),'r-');
    plot(xx,-hti,'k-','linewidth',2);
    plot(xx,z, 'k-','linewidth',2)
    plot(xx,TH(i,:).*z,'b--', 'linewidth',2); 
    ylim([ -dep 2]);
    drawnow
    system('sleep 0.1');
    waitbar(i/NT,h,'working...');
    if abort;break; end
end

delete(h)
