Up1 = [Un;Un(end,:)]; Up1 = .5*(Up1(1:end-1,:)+Up1(2:end,:));
Vp1 = [zeros(NX,1) Vn]; Vp1 = .5*(Vp1(:,2:end)+Vp1(:,1:end-1));

Up1(isnan(Up1))=0; Vp1(isnan(Vp1))=0;

INC = 2;
quiver(XXp(1:INC:end,1:INC:end),YYp(1:INC:end,1:INC:end), ...
       Up1(1:INC:end,1:INC:end), Vp1(1:INC:end,1:INC:end),'-b');
% 
% dy=yy(2)-yy(1);
% st = stream2(XXp',YYp',Up1',Vp1',-(B2-B1)+0.1+0*linspace(dy/2,A1,10), ...
%      linspace(dy/2,A1,10),[0.05 1000000]);
% streamline( st );
% 
% s1 = cell2mat(st(1,1));
% 
% xv = 0.5;
% [~,i]=min(abs(xv-s1(:,1)));
% aux = linspace(-(A2-A1),s1(i,2)-dy,15);
% 
% streamline( XXp',YYp',Up1',Vp1', xv+0*aux,aux, [0.001 80000]);