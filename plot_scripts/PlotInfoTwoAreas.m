
% plot Rho info measures
% two areas

syncthres = 0.90;%0.85; %0.4

cscat = coh; %peakcor;%
pscat = phasediff(:,:,1);
% pscat(phasediff_std > 0.1) = NaN;
pscat(cscat < syncthres) = NaN;

% cscat = peakcor(:,:,1);
% pscat = phasediff(:,:,1);%CVind_e(:,:,2); %phasediff(:,:,1);
% pscat(cscat < syncthres) = NaN;

clims = [-0.05,1.05]; %[0.04,0.16]; %[0,1];
rlims = [-0.2, 1.05];
ylims = [0, 1.5];

%% noise and signals

% % base
% figure;
% yscat3 = rho_base(:,:,3);
% yscat2 = rho_base(:,:,2);
% yscat = rho_base(:,:,1);
% subplot(122); scatter(pscat(:), yscat2(:), 30, sett.purple, 'filled');
% hold on; scatter(pscat(:), yscat(:), 30, sett.red, 'filled');
% hold on; scatter(pscat(:), yscat3(:), 30, sett.green, 'filled');
% xlabel('Phase difference'); ylabel('Correlation coefficient')
% xlim(xlims); ylim([-1,1]);
% subplot(121); scatter(cscat(:), yscat2(:), 30, sett.purple, 'filled');
% hold on; plot([syncthres, syncthres], [-1,1], ':k');
% scatter(cscat(:), yscat(:), 30, sett.green, 'filled');
% xlim(xlims);
% xlabel('Synchrony'); ylabel('Correlation coefficient')
% suptitle('Baseline firing rate');
% 
% 
% % power
% % rho_powermat = [rho_powermat, zeros(40,8,3)];
% figure; 
% yscat3 = rho_peak(:,:,3); % rho_power(:,:,3);
% yscat2 = rho_peak(:,:,2); % rho_power(:,:,2);
% yscat = rho_peak(:,:,1); % rho_power(:,:,1);
% subplot(122); scatter(pscat(:), yscat2(:), 30, sett.purple, 'filled');
% hold on; scatter(pscat(:), yscat(:), 30, sett.red, 'filled');
% hold on; scatter(pscat(:), yscat3(:), 30, sett.green, 'filled');
% xlabel('Phase difference'); ylabel('Correlation coefficient')
% xlim(xlims); ylim([-1,1]);
% subplot(121); scatter(cscat(:), yscat2(:), 30, sett.purple, 'filled');
% hold on; scatter(cscat(:), yscat(:), 30, sett.red, 'filled');
% hold on; scatter(cscat(:), yscat3(:), 30, sett.green, 'filled');
% hold on; plot([syncthres, syncthres], [-1,1], ':k');
% xlim(xlims);
% xlabel('Synchrony'); ylabel('Correlation coefficient')
% suptitle('LFP amplitude');
% 
% % freq
% figure; 
% yscat3 = rho_freq(:,:,3);
% yscat2 = rho_freq(:,:,2);
% yscat = rho_freq(:,:,1);
% subplot(122); scatter(pscat(:), yscat2(:), 30, sett.purple, 'filled');
% hold on; scatter(pscat(:), yscat(:), 30, sett.red, 'filled');
% hold on; scatter(pscat(:), yscat3(:), 30, sett.green, 'filled');
% xlabel('Phase difference'); ylabel('Correlation coefficient')
% xlim(xlims); ylim([-1,1]);
% subplot(121); scatter(cscat(:), yscat2(:), 30, sett.purple, 'filled');
% hold on; scatter(cscat(:), yscat(:), 30, sett.green, 'filled');
% hold on; plot([syncthres, syncthres], [-1,1], ':k');
% xlim(xlims);
% xlabel('Synchrony'); ylabel('Correlation coefficient')
% suptitle('Frequency');

%% other measures than sync and phase
% 
% Id = peakcor(:,:,1) > 0.4;
% 
% % CV
% figure
% subplot(121); hold on
% yscat2 = rho_power(:,:,2);
% yscat3 = rho_power(:,:,3);
% yscat = rho_power(:,:,1);
% zdata = CVe(:,:,1);
% scatter(zdata(Id), yscat2(Id), 30, sett.purple, 'filled');
% scatter(zdata(Id), yscat(Id), 30, sett.red, 'filled');
% scatter(zdata(Id), yscat3(Id), 30, sett.green, 'filled');
% xlim([0.04,0.16]); ylim([-1,1]);
% xlabel('Sync area 2');
% ylabel('Correlation coefficient')
% suptitle('LFP power')
% 
% % con
% subplot(122)
% zdata = repmat(0:0.02:0.3,[40,1]);
% scatter(zdata(Id), yscat(Id), 30, sett.red, 'filled');
% % xlim([0,1]); ylim([0.04, 0.16]);
% xlabel('Connection strength');
% ylabel('Correlation coefficient')
% 
% %

%% real cor

figure; 
yscat = rho_freq_12(:,:);

subplot(121); hold on; 
plot(clims, [0,0], 'k');
scatter(cscat(coh > syncthres), yscat(coh > syncthres), 30, sett.blue, 'filled');
scatter(cscat(coh < syncthres), yscat(coh < syncthres), 20, sett.blue);
plot([syncthres, syncthres], rlims, ':k');
xlim(clims); ylim(rlims);
xlabel('Coherence'); ylabel('Correlation coefficient')
[rho pval] = corr(cscat(:), yscat(:))

subplot(122); hold on; 
plot(clims, [0,0], 'k');
scatter(pscat(:), yscat(:), 30, sett.blue, 'filled');
xlabel('Phase difference'); ylabel('Correlation coefficient')
xlim(clims); ylim(rlims);

[rho pval] = circ_corrcl(pscat(~isnan(pscat))*pi*2, yscat(~isnan(pscat)))


suptitle('Oscillation frequency');


figure; 
yscat = rho_power_12(:,:);

subplot(121); hold on; 
plot(clims, [0,0], 'k');
scatter(cscat(coh > syncthres), yscat(coh > syncthres), 30, sett.green, 'filled');
scatter(cscat(coh < syncthres), yscat(coh < syncthres), 20, sett.green);
plot([syncthres, syncthres], rlims, ':k');
xlim(clims); ylim(rlims);
xlabel('Coherence'); ylabel('Correlation coefficient')
[rho pval] = corr(cscat(:), yscat(:))

subplot(122); hold on; 
plot(clims, [0,0], 'k');
scatter(pscat(:), yscat(:), 30, sett.green, 'filled');
xlabel('Phase difference'); ylabel('Correlation coefficient')
xlim(clims); ylim(rlims);
[rho pval] = circ_corrcl(pscat(~isnan(pscat))*2*pi, yscat(~isnan(pscat)))

suptitle('Power at oscillation frequency');

% weighted mean phase
dsel = ~isnan(pscat(:));
[sortyscat,Isort] = sort(yscat(dsel));
pscatn = pscat(dsel);
sortpscat = pscatn(Isort);
wmPhaseP = nanmean(sortpscat(end-9:end))
maxCorP = nanmean(sortyscat(end-9:end))


figure; 
yscat = rho_fper_12(:,:);

subplot(121); hold on; 
plot(clims, [0,0], 'k');
scatter(cscat(coh > syncthres), yscat(coh > syncthres), 30, sett.purple, 'filled');
scatter(cscat(coh < syncthres), yscat(coh < syncthres), 20, sett.purple);
plot([syncthres, syncthres], rlims, ':k');
xlim(clims); ylim(rlims);
xlabel('Coherence'); ylabel('Correlation coefficient')
[rho pval] = corr(cscat(:), yscat(:))

subplot(122); hold on; 
plot(clims, [0,0], 'k');
scatter(pscat(:), yscat(:), 30, sett.purple, 'filled');
xlabel('Phase difference'); ylabel('Correlation coefficient')
xlim(clims); ylim(rlims);
[rho pval] = circ_corrcl(pscat(~isnan(pscat))*pi*2, yscat(~isnan(pscat)))

suptitle('Firing rate');

% weighted mean phase
dsel = ~isnan(pscat(:));
[sortyscat,Isort] = sort(yscat(dsel));
pscatn = pscat(dsel);
sortpscat = pscatn(Isort);
wmPhaseF = mean(sortpscat(end-9:end))
maxCorF = mean(sortyscat(end-9:end))


%% peaks for individual runs

syncthres = 0.90;

cscat = coh;
pscat = phasediff(:,:,1);
pscat(cscat < syncthres) = NaN;
dsel = ~isnan(pscat(:));

nsel = 5;

for i = 1:5
    i
    
    %power
    yscat = rho_power_12LR(:,:,:,i);
    
    % weighted mean phase
    [sortyscat,Isort] = sort(yscat(dsel));
    pscatn = pscat(dsel);
    sortpscat = pscatn(Isort);
    wmPhaseP(i) = mean(sortpscat(end-(nsel-1):end));
    maxCorP(i) = mean(sortyscat(end-(nsel-1):end));
    
    
    % firing rate
    yscat = rho_fper_12LR(:,:,:,i);

    % weighted mean phase
    dsel = ~isnan(pscat(:));
    [sortyscat,Isort] = sort(yscat(dsel));
    pscatn = pscat(dsel);
    sortpscat = pscatn(Isort);
    wmPhaseF(i) = mean(sortpscat(end-(nsel-1):end));
    maxCorF(i) = mean(sortyscat(end-(nsel-1):end));
end

out = [maxCorP', wmPhaseP', maxCorF', wmPhaseF']

%% mutual info

figure; 
yscat = MI_freq(:,:);

subplot(121); hold on; 
plot(clims, [0,0], 'k');
scatter(cscat(coh > syncthres), yscat(coh > syncthres), 30, sett.blue, 'filled');
scatter(cscat(coh < syncthres), yscat(coh < syncthres), 20, sett.blue);
plot([syncthres, syncthres], rlims, ':k');
xlim(clims); 
ylim([0,1.05]);
xlabel('Coherence'); ylabel('Mutual information (bits)')
[rho pval] = corr(cscat(:), yscat(:))

subplot(122); hold on; 
plot(clims, [0,0], 'k');
scatter(pscat(:), yscat(:), 30, sett.blue, 'filled');
xlabel('Phase difference'); ylabel('Mutual information (bits)')
xlim(clims); 
ylim([0,1.05]);
[rho pval] = circ_corrcl(pscat(~isnan(pscat))*2*pi, yscat(~isnan(pscat)))

suptitle('Oscillation frequency');


figure; 
yscat = MI_power(:,:);

subplot(121); hold on; 
plot(clims, [0,0], 'k');
scatter(cscat(coh > syncthres), yscat(coh > syncthres), 30, sett.green, 'filled');
scatter(cscat(coh < syncthres), yscat(coh < syncthres), 20, sett.green);
plot([syncthres, syncthres], rlims, ':k');
xlim(clims); 
ylim([0, 0.3]);
xlabel('Coherence'); ylabel('Mutual information (bits)')
[rho pval] = corr(cscat(:), yscat(:))

subplot(122); hold on; 
plot(clims, [0,0], 'k');
scatter(pscat(:), yscat(:), 30, sett.green, 'filled');
xlabel('Phase difference'); ylabel('Mutual information (bits)')
xlim(clims); 
ylim([0, 0.3]);
[rho pval] = circ_corrcl(pscat(~isnan(pscat))*2*pi, yscat(~isnan(pscat)))

suptitle('Power at oscillation frequency');


figure; 
yscat = MI_fper(:,:);

subplot(121); hold on; 
plot(clims, [0,0], 'k');
scatter(cscat(coh > syncthres), yscat(coh > syncthres), 30, sett.purple, 'filled');
scatter(cscat(coh < syncthres), yscat(coh < syncthres), 20, sett.purple);
plot([syncthres, syncthres], rlims, ':k');
xlim(clims); 
ylim([0,0.15]);
xlabel('Coherence'); ylabel('Mutual information (bits)')
[rho pval] = corr(cscat(:), yscat(:))

subplot(122); hold on; 
plot(clims, [0,0], 'k');
scatter(pscat(:), yscat(:), 30, sett.purple, 'filled');
xlabel('Phase difference'); ylabel('Mutual information (bits)')
xlim(clims); 
ylim([0.1,0.25]);
[rho pval] = circ_corrcl(pscat(~isnan(pscat))*2*pi, yscat(~isnan(pscat)))

suptitle('Firing rate');


% corr and mi correlations

[rho pval] = corr(rho_fper_12(:), MI_fper(:))
[rho pval] = corr(rho_freq_12(:), MI_freq(:))
[rho pval] = corr(rho_power_12(:), MI_power(:))


%% for multiple runs in one plot...

    
    pos = [50,50,1500,500];
    
    syncthres = 0.90;
    
    % fper
    clims = [-0.05,1.05]; % coherence
    rlims = [-0.1, 0.55]; % rho
    
    pllims = [0,1];
    flims = [40,90];
    f2lims = [30,120];
    maxrn = 5;
    
    figure('Position', pos); hold on
    colors = repmat(sett.purple,[6,1])-0.1 + [([1:6]/9)', ([1:6]/9)', ([1:6]/9)'];
    for rn = 1:maxrn
        cscat = cohLR(:,:,1,rn);
        pscat = phasediffLR(:,:,1, rn);
        pscat(cscat < syncthres) = NaN;
        yscat = rho_fper_12LR(:,:,rn);
        plscat1 = PhaselockeLR(:,:,2,rn);
        plscat2 = Phaselocke12LR(:,:,1,rn);
        fscat1 = avgmorfreqLR(:,:,1,rn);
        fscat2 = avgmorfreqLR(:,:,2,rn);
        
        subplot(161); hold on;
        plot(clims, [0,0], 'k');
        scatter(cscat(coh > syncthres), yscat(coh > syncthres), 20, colors(rn,:), '*');
        scatter(cscat(coh < syncthres), yscat(coh < syncthres), 20, colors(rn,:));
        plot([syncthres, syncthres], rlims, ':k');

        subplot(162); hold on;
        plot(clims, [0,0], 'k');
        scatter(pscat(:), yscat(:), 20, colors(rn,:), '*');
        
        subplot(163); hold on;
        plot(clims, [0,0], 'k');
        scatter(plscat1(:), yscat(:), 10, colors(rn,:));
        
        subplot(164); hold on;
        plot(clims, [0,0], 'k');
        scatter(plscat2(:), yscat(:), 10, colors(rn,:));
        
        subplot(165); hold on;
        plot(clims, [0,0], 'k');
        scatter(fscat1(:), yscat(:), 10, colors(rn,:));
        
        subplot(166); hold on;
        plot(clims, [0,0], 'k');
        scatter(fscat2(:), yscat(:), 10, colors(rn,:));
        
    end
    subplot(161);
    xlim(clims); ylim(rlims);
    xlabel('Coherence'); ylabel('Correlation coefficient')
    
    subplot(162);
    xlabel('Phase difference');
    xlim(clims); ylim(rlims);
    
    subplot(163);
    xlabel('PPC_{2\rightarrow1}');
    xlim(pllims); ylim(rlims);
    dat = mean(PhaselockeLR(:,:,2,:),4);
    dat2 = mean(rho_fper_12LR,4);
    [rho, pval] = corr(dat(:), dat2(:))
    
    subplot(164);
    xlabel('PPC_{2\rightarrow2}');
    xlim(pllims); ylim(rlims);
    dat = mean(Phaselocke12LR(:,:,1,:),4);
   [rho, pval] = corr(dat(:), dat2(:))
    
    subplot(165);
    xlabel('freq_{1} (Hz)');
    xlim(flims); ylim(rlims);
    dat = mean(avgmorfreq(:,:,1,:),4);
    [rho, pval] = corr(dat(:), dat2(:))
    
    subplot(166);
    xlabel('freq_{2} (Hz)');
    xlim(flims); ylim(rlims);
    dat = mean(avgmorfreq(:,:,2,:),4);
    [rho, pval] = corr(dat(:), dat2(:))

    
    % freq
    rlims = [-0.5, 1.05]; % rho
figure('Position', pos); hold on;
    colors = parula(10);
    for rn = 1:maxrn
        cscat = cohLR(:,:,1,rn);
        pscat = phasediffLR(:,:,1, rn);
        pscat(cscat < syncthres) = NaN;
        yscat = rho_freq_12LR(:,:,rn);
        plscat1 = PhaselockeLR(:,:,2,rn);
        plscat2 = Phaselocke12LR(:,:,1,rn);
        fscat1 = avgmorfreqLR(:,:,1,rn);
        fscat2 = avgmorfreqLR(:,:,2,rn);
        
        subplot(161); hold on;
        plot(clims, [0,0], 'k');
        scatter(cscat(coh > syncthres), yscat(coh > syncthres), 20, colors(rn,:), '*');
        scatter(cscat(coh < syncthres), yscat(coh < syncthres), 20, colors(rn,:));
        plot([syncthres, syncthres], rlims, ':k');

        subplot(162); hold on;
        plot(clims, [0,0], 'k');
        scatter(pscat(:), yscat(:), 20, colors(rn,:), '*');
                
        subplot(163); hold on;
        plot(clims, [0,0], 'k');
        scatter(plscat1(:), yscat(:), 10, colors(rn,:));
        
        subplot(164); hold on;
        plot(clims, [0,0], 'k');
        scatter(plscat2(:), yscat(:), 10, colors(rn,:));
        
        subplot(165); hold on;
        plot(clims, [0,0], 'k');
        scatter(fscat1(:), yscat(:), 10, colors(rn,:));
        
        subplot(166); hold on;
        plot(clims, [0,0], 'k');
        scatter(fscat2(:), yscat(:), 10, colors(rn,:));
        
    end
    subplot(161);
    xlim(clims); ylim(rlims);
    xlabel('Coherence'); ylabel('Correlation coefficient')
    
    subplot(162);
    xlabel('Phase difference');
    xlim(clims); ylim(rlims);
    
    subplot(163);
    xlabel('PPC_{2\rightarrow1}');
    xlim(pllims); ylim(rlims);
    dat = mean(PhaselockeLR(:,:,2,:),4);
    dat2 = mean(rho_freq_12LR,4);
    [rho, pval] = corr(dat(:), dat2(:))
    
    subplot(164);
    xlabel('PPC_{2\rightarrow2}');
    xlim(pllims); ylim(rlims);
    dat = mean(Phaselocke12LR(:,:,1,:),4);
   [rho, pval] = corr(dat(:), dat2(:))
    
    subplot(165);
    xlabel('freq_{1} (Hz)');
    xlim(flims); ylim(rlims);
    
    subplot(166);
    xlabel('freq_{2} (Hz)');
    xlim(flims); ylim(rlims);
    
    
    % power
    figure('Position', pos); hold on
    rlims = [-0.4, 1.05]; % rho
        colors = summer(8);
    for rn = 1:maxrn
        cscat = cohLR(:,:,1,rn);
        pscat = phasediffLR(:,:,1, rn);
        pscat(cscat < syncthres) = NaN;
        yscat = rho_power_12LR(:,:,rn);
        plscat1 = PhaselockeLR(:,:,2,rn);
        plscat2 = Phaselocke12LR(:,:,1,rn);
        fscat1 = avgmorfreqLR(:,:,1,rn); %reLR(:,:,2,rn); %mean(reLR(:,:,2,:),4);
        fscat2 = avgmorfreqLR(:,:,2,rn);
        
        subplot(161); hold on;
        plot(clims, [0,0], 'k');
        scatter(cscat(coh > syncthres), yscat(coh > syncthres), 20, colors(rn,:), '*');
        scatter(cscat(coh < syncthres), yscat(coh < syncthres), 20, colors(rn,:));
        plot([syncthres, syncthres], rlims, ':k');

        subplot(162); hold on;
        plot(clims, [0,0], 'k');
        scatter(pscat(:), yscat(:), 20, colors(rn,:), '*');
        
        subplot(163); hold on;
        plot(clims, [0,0], 'k');
        scatter(plscat1(:), yscat(:), 10, colors(rn,:));
        
        subplot(164); hold on;
        plot(clims, [0,0], 'k');
        scatter(plscat2(:), yscat(:), 10, colors(rn,:));
        
        subplot(165); hold on;
        plot(clims, [0,0], 'k');
        scatter(fscat1(:), yscat(:), 10, colors(rn,:));
        
        subplot(166); hold on;
        plot(clims, [0,0], 'k');
        scatter(fscat2(:), yscat(:), 10, colors(rn,:));
        
    end
    subplot(161);
    xlim(clims); ylim(rlims);
    xlabel('Coherence'); ylabel('Correlation coefficient')
    
    subplot(162);
    xlabel('Phase difference');
    xlim(clims); ylim(rlims);
    
    subplot(163);
    xlabel('PPC_{2\rightarrow1}');
    xlim(pllims); ylim(rlims);
    dat = mean(PhaselockeLR(:,:,2,:),4);
    dat2 = mean(rho_power_12LR,4);
    [rho, pval] = corr(dat(:), dat2(:))
    
    subplot(164);
    xlabel('PPC_{2\rightarrow2}');
    xlim(pllims); ylim(rlims);
    dat = mean(Phaselocke12LR(:,:,1,:),4);
   [rho, pval] = corr(dat(:), dat2(:))
    
    subplot(165);
    xlabel('Firing rate_{2} (Hz)');
    xlim(flims); ylim(rlims);
%     dat = mean(reLR(:,:,2,:),4);
%     figure; scatter(dat(:), dat2(:))
%    [rho, pval] = corr(dat(:), dat2(:))
    
    subplot(166);
    xlabel('freq_{2} (Hz)');
    xlim(flims); ylim(rlims);
    
    
    %%
    
    % power
    figure('Position', pos); hold on
    rlims = [-0.4, 1.05]; % rho
        colors = summer(8);
    for rn = 1:maxrn
        cscat = cohLR(:,:,1,rn);
        pscat = phasediffLR(:,:,1, rn);
        pscat(cscat < syncthres) = NaN;
        yscat = MI_power_12LR(:,:,rn);
        plscat1 = PhaselockeLR(:,:,2,rn);
        plscat2 = Phaselocke12LR(:,:,1,rn);
        fscat1 = mean(reLR(:,:,2,:),4);%reLR(:,:,2,rn);
        fscat2 = avgmorfreqLR(:,:,2,rn);
        
        subplot(161); hold on;
        plot(clims, [0,0], 'k');
        scatter(cscat(coh > syncthres), yscat(coh > syncthres), 20, colors(rn,:), '*');
        scatter(cscat(coh < syncthres), yscat(coh < syncthres), 20, colors(rn,:));
        plot([syncthres, syncthres], rlims, ':k');

        subplot(162); hold on;
        plot(clims, [0,0], 'k');
        scatter(pscat(:), yscat(:), 20, colors(rn,:), '*');
        
        subplot(163); hold on;
        plot(clims, [0,0], 'k');
        scatter(plscat1(:), yscat(:), 10, colors(rn,:));
        
        subplot(164); hold on;
        plot(clims, [0,0], 'k');
        scatter(plscat2(:), yscat(:), 10, colors(rn,:));
        
        subplot(165); hold on;
        plot(clims, [0,0], 'k');
        scatter(fscat1(:), yscat(:), 10, colors(rn,:));
        
        subplot(166); hold on;
        plot(clims, [0,0], 'k');
        scatter(fscat2(:), yscat(:), 10, colors(rn,:));
        
    end
    subplot(161);
    xlim(clims); ylim(rlims);
    xlabel('Coherence'); ylabel('Correlation coefficient')
    
    subplot(162);
    xlabel('Phase difference');
    xlim(clims); ylim(rlims);
    
    subplot(163);
    xlabel('PPC_{2\rightarrow1}');
    xlim(pllims); ylim(rlims);
    
    subplot(164);
    xlabel('PPC_{2\rightarrow2}');
    xlim(pllims); ylim(rlims);
    
    subplot(165);
    xlabel('Firing rate_{2} (Hz)');
    xlim(f2lims); ylim(rlims);
    
    subplot(166);
    xlabel('freq_{2} (Hz)');
    xlim(flims); ylim(rlims);
    
    
%% shuffled correlations

% time shuffled
figure; 
yscat = rho_freq_12t(:,:);

subplot(121); hold on; 
plot(clims, [0,0], 'k');
scatter(cscat(coh > syncthres), yscat(coh > syncthres), 30, sett.blue, 'filled');
scatter(cscat(coh < syncthres), yscat(coh < syncthres), 20, sett.blue);
plot([syncthres, syncthres], rlims, ':k');
xlim(clims); ylim(rlims);
xlabel('Coherence'); ylabel('Correlation coefficient')
[rho pval] = corr(cscat(:), yscat(:))

subplot(122); hold on; 
plot(clims, [0,0], 'k');
scatter(pscat(:), yscat(:), 30, sett.blue, 'filled');
xlabel('Phase difference'); ylabel('Correlation coefficient')
xlim(clims); ylim(rlims);

[rho pval] = circ_corrcl(pscat(~isnan(pscat))*pi*2, yscat(~isnan(pscat)))


suptitle('Oscillation frequency');


figure; 
yscat = rho_power_12t(:,:);

subplot(121); hold on; 
plot(clims, [0,0], 'k');
scatter(cscat(coh > syncthres), yscat(coh > syncthres), 30, sett.green, 'filled');
scatter(cscat(coh < syncthres), yscat(coh < syncthres), 20, sett.green);
plot([syncthres, syncthres], rlims, ':k');
xlim(clims); ylim(rlims);
xlabel('Coherence'); ylabel('Correlation coefficient')
[rho pval] = corr(cscat(:), yscat(:))

subplot(122); hold on; 
plot(clims, [0,0], 'k');
scatter(pscat(:), yscat(:), 30, sett.green, 'filled');
xlabel('Phase difference'); ylabel('Correlation coefficient')
xlim(clims); ylim(rlims);
[rho pval] = circ_corrcl(pscat(~isnan(pscat))*2*pi, yscat(~isnan(pscat)))

suptitle('Power at oscillation frequency');

% weighted mean phase
dsel = ~isnan(pscat(:));
[sortyscat,Isort] = sort(yscat(dsel));
pscatn = pscat(dsel);
sortpscat = pscatn(Isort);
wmPhaseP = nanmean(sortpscat(end-9:end))
maxCorP = nanmean(sortyscat(end-9:end))


figure; 
yscat = rho_fper_12t(:,:);

subplot(121); hold on; 
plot(clims, [0,0], 'k');
scatter(cscat(coh > syncthres), yscat(coh > syncthres), 30, sett.purple, 'filled');
scatter(cscat(coh < syncthres), yscat(coh < syncthres), 20, sett.purple);
plot([syncthres, syncthres], rlims, ':k');
xlim(clims); ylim(rlims);
xlabel('Coherence'); ylabel('Correlation coefficient')
[rho pval] = corr(cscat(:), yscat(:))

subplot(122); hold on; 
plot(clims, [0,0], 'k');
scatter(pscat(:), yscat(:), 30, sett.purple, 'filled');
xlabel('Phase difference'); ylabel('Correlation coefficient')
xlim(clims); ylim(rlims);
[rho pval] = circ_corrcl(pscat(~isnan(pscat))*pi*2, yscat(~isnan(pscat)))

suptitle('Firing rate');

% weighted mean phase
dsel = ~isnan(pscat(:));
[sortyscat,Isort] = sort(yscat(dsel));
pscatn = pscat(dsel);
sortpscat = pscatn(Isort);
wmPhaseF = mean(sortpscat(end-9:end))
maxCorF = mean(sortyscat(end-9:end))
   


% period shuffled
figure; 
yscat = rho_freq_12p(:,:);

subplot(121); hold on; 
plot(clims, [0,0], 'k');
scatter(cscat(coh > syncthres), yscat(coh > syncthres), 30, sett.blue, 'filled');
scatter(cscat(coh < syncthres), yscat(coh < syncthres), 20, sett.blue);
plot([syncthres, syncthres], rlims, ':k');
xlim(clims); ylim(rlims);
xlabel('Coherence'); ylabel('Correlation coefficient')
[rho pval] = corr(cscat(:), yscat(:))

subplot(122); hold on; 
plot(clims, [0,0], 'k');
scatter(pscat(:), yscat(:), 30, sett.blue, 'filled');
xlabel('Phase difference'); ylabel('Correlation coefficient')
xlim(clims); ylim(rlims);

[rho pval] = circ_corrcl(pscat(~isnan(pscat))*pi*2, yscat(~isnan(pscat)))


suptitle('Oscillation frequency');


figure; 
yscat = rho_power_12p(:,:);

subplot(121); hold on; 
plot(clims, [0,0], 'k');
scatter(cscat(coh > syncthres), yscat(coh > syncthres), 30, sett.green, 'filled');
scatter(cscat(coh < syncthres), yscat(coh < syncthres), 20, sett.green);
plot([syncthres, syncthres], rlims, ':k');
xlim(clims); ylim(rlims);
xlabel('Coherence'); ylabel('Correlation coefficient')
[rho pval] = corr(cscat(:), yscat(:))

subplot(122); hold on; 
plot(clims, [0,0], 'k');
scatter(pscat(:), yscat(:), 30, sett.green, 'filled');
xlabel('Phase difference'); ylabel('Correlation coefficient')
xlim(clims); ylim(rlims);
[rho pval] = circ_corrcl(pscat(~isnan(pscat))*2*pi, yscat(~isnan(pscat)))

suptitle('Power at oscillation frequency');

% weighted mean phase
dsel = ~isnan(pscat(:));
[sortyscat,Isort] = sort(yscat(dsel));
pscatn = pscat(dsel);
sortpscat = pscatn(Isort);
wmPhaseP = nanmean(sortpscat(end-9:end))
maxCorP = nanmean(sortyscat(end-9:end))


figure; 
yscat = rho_fper_12p(:,:);

subplot(121); hold on; 
plot(clims, [0,0], 'k');
scatter(cscat(coh > syncthres), yscat(coh > syncthres), 30, sett.purple, 'filled');
scatter(cscat(coh < syncthres), yscat(coh < syncthres), 20, sett.purple);
plot([syncthres, syncthres], rlims, ':k');
xlim(clims); ylim(rlims);
xlabel('Coherence'); ylabel('Correlation coefficient')
[rho pval] = corr(cscat(:), yscat(:))

subplot(122); hold on; 
plot(clims, [0,0], 'k');
scatter(pscat(:), yscat(:), 30, sett.purple, 'filled');
xlabel('Phase difference'); ylabel('Correlation coefficient')
xlim(clims); ylim(rlims);
[rho pval] = circ_corrcl(pscat(~isnan(pscat))*pi*2, yscat(~isnan(pscat)))

suptitle('Firing rate');

% weighted mean phase
dsel = ~isnan(pscat(:));
[sortyscat,Isort] = sort(yscat(dsel));
pscatn = pscat(dsel);
sortpscat = pscatn(Isort);
wmPhaseF = mean(sortpscat(end-9:end))
maxCorF = mean(sortyscat(end-9:end))