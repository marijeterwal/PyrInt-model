
%% set thresholds for visualization

freq1 = 10;
freq2 = 40;
freq3 = 60;
rratio = 0.75;

ple0 = 0.05;%0.025;
ple1 = 0.25;%0.075;
ple2 = 0.50;%0.115;
pli0 = 0.05; %0.18

re0 = 0.1;
ri0 = 0.1;

plotPle(meanRe < re0) = NaN;
plotPli(meanRi < re0) = NaN;
plotPle(plotPle <= ple0 | plotPli <= pli0) = NaN;
plotPli(plotPle <= ple0 | plotPli <= pli0) = NaN;

plotPle(isnan(plotPle)) = 0;
plotPli(isnan(plotPli)) = 0;

plotFreqe(plotPle <= ple0) = NaN;
stdFreq(plotPle <= ple0) = NaN;
plotFreqi(plotPli <= pli0) = NaN;
stdFreqi(plotPli <= pli0) = NaN;

%% plot image of effect of drive - one area

plotMat = zeros(l1steps, l2steps);

plotMat(plotPli <= pli0 | plotPle <= ple0) = 0; % async
plotMat(plotPle <= ple1 & plotFreqe > freq1 & plotFreqe <= freq2) = 1; %lowsyn lowgam
plotMat(plotPle <= ple2 & plotPle > ple1 & plotFreqe > freq1 & plotFreqe <= freq2) = 2; %medsyn lowgam
plotMat(plotPle > ple2 & plotFreqe > freq1 & plotFreqe <= freq2) = 3; %highsyn lowgam
plotMat(plotPle <= ple1 & plotFreqe > freq2 & plotFreqe <= freq3) = 4; %lowsyn medgam
plotMat(plotPle <= ple2 & plotPle > ple1 & plotFreqe > freq2 & plotFreqe <= freq3) = 5; %medsyn medgam
plotMat(plotPle > ple2 & plotFreqe > freq2 & plotFreqe <= freq3) = 6; %highsyn medgam
plotMat(plotPle <= ple1 & plotFreqe > freq3) = 7; %lowsyn medgam
plotMat(plotPle <= ple2 & plotPle > ple1 & plotFreqe > freq3) = 8; %medsyn medgam
plotMat(plotPle > ple2 & plotFreqe > freq3) = 9; %highsyn medgam

plotMat(plotPli > pli0 & plotPle < ple0 & meanRe < re0) = 11; % ING no pyr

figure ('Name','Effect of drive - summary', 'Position',[100 100 500 400]);
map = [ [0.8,0.8,0.8]%[0 0 143]/255; ...
        [179 251 255]/255; [114 232 226]/255; [1 190 199]/255; ...
        [255 233 171]/255; [255 220 121]/255; [255 192 13]/255; ...
        [255 173 171]/255; [254 129 126]/255; [254 80 76]/255; ...
        [0.6,0.6,0.6]; [0.4,0.4,0.4]]; %[40 40 203]/255];%

xax = sett.iIu:sett.Istep:sett.iIu_max;
yax = sett.eIu:sett.Istep:sett.eIu_max;

xl = [sett.iIu(1), sett.iIu_max(1)];
yl =  [sett.eIu(1), sett.eIu_max(1)];

imagesc(xax, yax, plotMat');
axis square;
caxis([0 11]);
axis image;
xlim([0,3])
set(gca,'YDir','normal');
colormap(map);
colorbar;
% title('');
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');

%% contours of phase locking
% figure; 
% subplot(121);
% contour(xax, yax, plotPle',0.025, 'color', sett.red, 'linewidth', 2)
% title('Pyramidal cells')
% xlabel('I_{i, inj} (\muA/cm^{2})');
% ylabel('I_{e, inj} (\muA/cm^{2})');
% subplot(122);
% contour(xax, yax, plotPli',0.15, 'color', sett.blue, 'linewidth', 2)
% title('Interneurons')
% xlabel('I_{i, inj} (\muA/cm^{2})');

%% freq and sync and their stds - when loading multiple runs with different seeds

% SI plot
figure;

subplot(231);
imagesc(xax, yax, plotFreqe)
axis xy
axis image
colormap(parula)
colorbar; caxis([20,100])
xlim(xl)
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');
title('Frequency - pyramidal cells')

subplot(232);
imagesc(xax, yax, plotPle)
axis xy
axis image
colormap(parula)
colorbar; caxis([0,1])
xlim(xl)
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');
title('PPC - pyramidal cells')

subplot(233);
imagesc(xax, yax, meanRe)
axis xy
axis image
colormap(parula)
colorbar; caxis([0,100])
xlim(xl)
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');
title('Firing rate - pyramidal cells')

subplot(234);
imagesc(xax, yax, plotFreqi)
axis xy
axis image
colormap(parula)
xlim(xl)
colorbar; caxis([20,100])
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');
title('Frequency - interneurons')

subplot(235);
imagesc(xax, yax, plotPli)
axis xy
axis image
colormap(parula)
xlim(xl)
colorbar; caxis([0,1])
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');
title('PPC - interneurons')

subplot(236);
imagesc(xax, yax, meanRi)
axis xy
axis image
colormap(parula)
colorbar; caxis([0,100])
xlim(xl)
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');
title('Firing rate - interneurons')


%frequency
figure;

meanFreq = plotFreqe;
meanFreqi = plotFreqi;

subplot(221);
imagesc(xax, yax, meanFreq)
axis xy
axis image
colorbar; caxis([20,150])
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');
title('Mean')

if exist('stdFreq','var')
subplot(222);
imagesc(xax, yax, stdFreq)
axis xy
axis image
colorbar; caxis([0,15])
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');
title('Standard deviation')
end

subplot(223);
imagesc(xax, yax, meanFreqi)
axis xy
axis image
colorbar; caxis([20,150])
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');
title('Mean')

if exist('stdFreqi','var')
subplot(224);
imagesc(xax, yax, stdFreqi)
axis xy
axis image
colorbar; caxis([0,15])
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');
title('Standard deviation')
end

% phase locking

figure;

subplot(221);
imagesc(xax, yax, plotPle)
axis xy
axis image
colorbar; caxis([0,0.20]);
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');
title('Mean')

if exist('Ple_std','var')
subplot(222);
imagesc(xax, yax, Ple_std)
axis xy
axis image
colorbar; caxis([0,0.1]);
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');
title('Standard deviation')
end

subplot(223);
imagesc(xax, yax, plotPli)
axis xy
axis image
colorbar; caxis([0,0.80]);
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');
title('Mean')

if exist('Pli_std','var')
subplot(224);
imagesc(xax, yax, Pli_std)
axis xy
axis image
colorbar; caxis([0,0.2]);
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');
title('Standard deviation')
end

% firing rate

figure;

subplot(221);
imagesc(xax, yax, meanRe)
axis xy
axis image
colorbar; caxis([0,150]);
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');
title('Mean')

if exist('stdRe','var')
subplot(222);
imagesc(xax, yax, stdRe)
axis xy
axis image
colorbar; caxis([0,5]);
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');
title('Standard deviation')
end

subplot(223);
imagesc(xax, yax, meanRi)
axis xy
axis image
colorbar; caxis([0,150]);
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');
title('Mean')

if exist('stdRi','var')
subplot(224);
imagesc(xax, yax, stdRi)
axis xy
axis image
colorbar; caxis([0,5]);
xlabel('I_{i, inj} (\muA/cm^{2})');
ylabel('I_{e, inj} (\muA/cm^{2})');
title('Standard deviation')
end

%% Single dimensions
figure; 
subplot(211);
plot(avgmorfreq(1:end), 'color', sett.purple, 'linewidth', 2);
ylim([30, 70]); xlim([1,26]);
subplot(212);
plot(Phaselocke(1:end), 'color', sett.green, 'linewidth', 2);
ylim([0.0, 0.6]);xlim([1,26]);