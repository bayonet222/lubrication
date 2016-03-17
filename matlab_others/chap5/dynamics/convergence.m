%% ANTES EJECUTAR main_re_ea_analitic.m

N = [64 128 256 512 1024];

errp = zeros(length(N),1);
errt = zeros(length(N),1);
errT = zeros(length(N),1);

for in = 1:length(N)
    
    dnu = zeros(N(in),1); dex = zeros(N(in),1);
    thetaex = ones(N(in),1);
    
    dx = 1/N(in);
    dataea = load(strcat(num2str(N(in),'%d'),'.dat'));
    dataea(:,1) = dataea(:,1)*10-0.5;
    dataea(:,2) = dataea(:,2)*10;
    
    I = find(dataea(:,1) >= 0.0 & dataea(:,1) < 1.0);
    dataea = dataea(I,:);
    peaexact = interp1([0 x1n beta(it) xx(end)],[0 pressea(it) 0 0],dataea(:,1));
    
    paux = find(dataea(:,5)==-1,1,'first'); %% LEFT SIDE OF THE POCKET
    
    paux1 = find(dataea(:,5)==-1,1,'last'); %% RIGHT SIDE OF THE POCKET
    [~,paux2] = min(abs(dataea(:,1)-beta(it))); %% RIGHT SIDE OF THE POCKET
    
    for j = 2:N-1
        dex(j) = (peaexact(j+1)-peaexact(j-1))/(2*dx);
        dnu(j) = (dataea(j+1,2)-dataea(j-1,2))/(2*dx);
    end
    dnu(paux) = dex(paux);
    dnu(paux-1) = dex(paux-1);
    
    thetaex(paux2:paux1) = 0.5; %% h1/h2
    
    errt(in) = sqrt(dx)*(norm(thetaex-dataea(:,3),2));
    errp(in) = sqrt(dx)*(norm(peaexact-dataea(:,2),2) + norm(dex-dnu,2));
end

figure(1)
plot(log2(errt),'-o')
set(gcf,'Color','w'); set(gca,'FontSize',14)

figure(2)
plot(log2(errp),'-o')
set(gcf,'Color','w'); set(gca,'FontSize',14)
