%% caso malla geometria del proyecto
II = reshape(1:NX*NY, NX, NY); % indices de todos los volumenes
IW = II(1,2:end-1)'; IE = II(end,2:end-1)'; 
IN = II(2:end-1,end); IS = II(2:end-1,1); %indices, oeste, este, norte, sur
INE = II(end,end); INW = II(1,end); ISE = II(end,1); ISW = II(1,1);

II([IW;IE;IN;IS;INE;INW;ISE;ISW]) = -1; % marco los nodos no internos

II = II(:); II(II==-1) = []; %quito los nodos no internos

%%
IP = sort([II;IN;IS;IW;IE;INW;INE;ISW;ISE]);%P(INW)=0
IPU = sort([II;IN;IS;IE;INE;ISE]);
IPV = sort([II;IS;ISW;IW;IE;ISE]);