
figure(num_fig);clf

subplot(2,1,1);
% 
pcolor([XXp-dx/2;XXp(end,:)+dx/2],[YYp;YYp(end,:)],[Un;Un(end,:)]); 

shading 'interp';colorbar;

subplot(2,1,2);hold on

VWALL = [0 0.5*(Vn(1,:)+Vn(end,:))];
pcolor([zeros(1,NY+1); XXp(:,1) XXp; zeros(1,NY+1)+L], [linspace(0,H,NY+1);YYp(:,1)-dx/2 YYp+dx/2; linspace(0,H,NY+1)], [VWALL;zeros(NX,1) Vn;VWALL]);
xlim([xx(1) xx(end)]);ylim([yy(1) yy(end)]); title('V');
shading 'interp';colorbar
graficas_campo_vel

% hold on;
% graficas_campo_vel
% colorbar; shading 'interp';

% subplot(3,1,3);
% hold on
% pcolor(XXp,YYp, Pn); 
% contourf(XXp,YYp, Pn,50);
% title(strcat( 'Re = ', num2str(Re,2),' tn = ', num2str(niter*dt,2),...
%        ' resU = ', num2str(resU,3), ' resP = ', num2str(resP,3),...
%        ' resG = ', num2str(resG,3) ));
%  shading 'interp';
% colorbar;   

drawnow

