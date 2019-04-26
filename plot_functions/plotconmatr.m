function plotconmatr(sett,conii, conie, conee, conei, saveind, saveloc)

reg = sett.Nregions;
Ni = sett.Ni;
Ne = sett.Ne;
N = Ni + Ne;

%% construct connection matrix
con = zeros(N*reg, N*reg);

if reg == 1; 
    ncomb = 1;
    comb = [1,1];
else
    ncomb = 2*nchoosek(reg,2);
    comb = [combnk(1:reg,2); fliplr(combnk(1:reg,2))];
end

for r = 1:reg
    con((r-1)*N+1:(r-1)*N+Ni,(r-1)*N+1:(r-1)*N+Ni)  = 1*    conii((r-1)*Ni+1:r*Ni,(r-1)*Ni+1:r*Ni);
    con((r-1)*N+Ni+1:r*N,(r-1)*N+1:(r-1)*N+Ni)      = 0.5*  conie((r-1)*Ne+1:r*Ne,(r-1)*Ni+1:r*Ni);
    con((r-1)*N+1:(r-1)*N+Ni, (r-1)*N+Ni+1:r*N)     = -0.5* conei((r-1)*Ni+1:r*Ni,(r-1)*Ne+1:r*Ne);
    con((r-1)*N+Ni+1:r*N, (r-1)*N+Ni+1:r*N)         = -1*   conee((r-1)*Ne+1:r*Ne,(r-1)*Ne+1:r*Ne);
end
for nc = 1:ncomb;
    reg1 = comb(nc,1); % from
    reg2 = comb(nc,2); % to
    
    if sett.g_syn_ee_r(reg2,reg1) ~= 0;
        con((reg2-1)*N+1:(reg2-1)*N+Ni, (reg1-1)*N+Ni+1:reg1*N) = -0.5*conei((reg2-1)*Ni+1:reg2*Ni,(reg1-1)*Ne+1:reg1*Ne);
    end
    if sett.g_syn_ei_r(reg2,reg1) ~= 0;
        con((reg2-1)*N+Ni+1:reg2*N, (reg1-1)*N+Ni+1:reg1*N) = -1*conee((reg2-1)*Ne+1:reg2*Ne,(reg1-1)*Ne+1:reg1*Ne);
    end
    
end

%% plot!
figure('Name','Connection matrix', 'Position',[100 200 800 600]); % temp
imagesc(1:1:N*reg, 1:1:N*reg, con);

xlim ([0.5 N*reg+0.5]); ylim ([0.5 N*reg+0.5]);
for i = 1:reg-1
    hold on;
    x = [0.5 N*reg+0.5]; y = [N*i+0.5 N*i+0.5]; plot(x,y, 'color', [0.1 0.1 0.1], 'LineWidth', 1.5);
    x = [N*i+0.5 N*i+0.5]; y = [0.5 N*reg+0.5]; plot(x,y, 'color', [0.1 0.1 0.1], 'LineWidth', 1.5);
end
if sett.subregion > 0;
    % Nisub and Nesub
    if sett.subregion_type == 0; Nisub = 0;
    elseif sett.subregion_type == 1 || sett.subregion_type == 2;
        Nisub = round(Ni*sett.subregion_numneurons/100);
    end
    if sett.subregion_type == 1; Nesub = 0;
    elseif sett.subregion_type == 0 || sett.subregion_type == 2
        Nesub = round(Ne*sett.subregion_numneurons/100);
    end
    % Lines
    for k = 1:sett.subregion;
        hold on;
        mid = N*(k-1)+ Ni;
        grijs = [0.4 0.4 0.4];
        x = [0.5 N*reg+0.5]; y = [mid-Nisub+0.5  mid-Nisub+0.5]; plot(x,y, 'color', grijs, 'LineWidth', 1);
        x = [0.5 N*reg+0.5]; y = [mid+Nesub+0.5 mid+Nesub+0.5]; plot(x,y, 'color', grijs, 'LineWidth', 1);
        x = [mid-Nisub+0.5 mid-Nisub+0.5]; y = [0.5 N*reg+0.5]; plot(x,y, 'color', grijs, 'LineWidth', 1);
        x = [mid+Nesub+0.5 mid+Nesub+0.5]; y = [0.5 N*reg+0.5]; plot(x,y, 'color', grijs, 'LineWidth', 1);
    end
end
axis image;
map = [[140 20 10]/256; sett.red;  [1 1 1];  sett.blue; [50 80 150]/256];
colormap(map);
title ('Connection matrix');
xlabel ('From'); ylabel ('To');
colorbar('YTick', [-0.8 -0.4 0 0.4 0.8],...
    'YTickLabel',...
    { 'E to E','E to I','No con', 'I to E', 'I to I'});
hold off;

if saveind == true
    saveas(gcf, [saveloc,'connectionmatrix'],'fig');
end
end