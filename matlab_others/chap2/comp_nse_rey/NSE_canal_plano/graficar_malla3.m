figure(777)
hold on
plot([XXp-dx/2;XXp(end,:)+dx/2],[YYp;YYp(end,:)],'xr',[XXp-dx/2;XXp(end,:)+dx/2]',[YYp;YYp(end,:)]','xr');
plot([XXp(:,1) XXp],[YYp(:,1)-dx/2 YYp+dx/2],'>b', [XXp(:,1) XXp]',[YYp(:,1)-dx/2 YYp+dx/2]','>b')
plot(XXp,YYp,'ok',XXp',YYp','ok');

[YY,XX]=meshgrid(linspace(0,H,NY+1),linspace(-L,L,NX+1));
plot(XX,YY,'-k',XX',YY','-k');

plot(XXp(II),YYp(II),'*r','linewidth',4);
plot(XXp(IS),YYp(IS),'*g','linewidth',4);
plot(XXp(IW),YYp(IW),'*y','linewidth',4);

plot(XXp(ISW),YYp(ISW),'*k','linewidth',4);
plot(XXp(IN),YYp(IN),'*m','linewidth',4);
plot(XXp(IEw),YYp(IEw),'*b','linewidth',4); plot(XXp(INEw),YYp(INEw),'xk','linewidth',2);
plot(XXp(IE),YYp(IE),'+m','linewidth',2);
% plot(XXp(IPU),YYp(IPU),'^y','linewidth',2);