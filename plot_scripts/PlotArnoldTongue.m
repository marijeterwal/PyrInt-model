
% plot image of effect of drive - one area

freq1 = 10;
freq2 = 40;
freq3 = 60;
rratio = 0.75;

ple0 = 0.025;
ple1 = 0.075;
ple2 = 0.115;

pli0 = 0.15;

xax = sett.Iper_freq:sett.perstep:sett.Iper_freq_max;
yax = sett.Iper_amp:sett.perampstep:sett.Iper_amp_max;

xlab = 'Frequency of drive (Hz)';
ylab = 'Amplitude of drive (a.u.)';

%% freq and sync and their stds

% meanFreq(plotPle < 0.01) = NaN;
% stdFreq(plotPle < 0.01) = NaN;
% meanFreqi(plotPli < 0.01) = NaN;
% stdFreqi(plotPli < 0.01) = NaN;

%frequency
figure;

fdat = repmat(xax',[1,length(yax)]);
favg = mean(meanFreq(:,1));

subplot(221);
% plotFr = 1-2*abs(0.5 - meanFreq.^1 ./(meanFreq.^1+repmat(xax',[1,length(yax)]).^1));
% plotFr = 1-sqrt(abs(0.5 - meanFreq.^2 ./(meanFreq.^2+repmat(xax',[1,length(yax)]).^2)));
% plotFr = 1-abs(meanFreq-fdat)./(abs(repmat(meanFreq(1,:),[length(xax),1])-fdat));
plotFr = 1-abs(meanFreq-fdat)./(xax(end)-favg);% (2*abs(repmat(35.1,[length(xax),length(yax)])-fdat));
binFr = zeros(size(plotFr));
% binFr(plotFr > 0.95 & plotFr < 1.05) = 1;
binFr(plotFr > 0.90) = 1;

imagesc(xax, yax, plotFr')
axis xy
axis square
colormap(parula)
hold on
contour(xax, yax' ,binFr',[0.5,0.5], 'w')
colorbar; caxis([0,1])
xlabel(xlab);
ylabel(ylab);
title('Mean')

subplot(222);
% plotFr = stdFreq./repmat(xax',[1,length(yax)]);
plotFr = 1-2*abs(0.5 - stdFreq.^2 ./(stdFreq.^2+repmat(xax',[1,length(yax)]).^2));
imagesc(xax, yax, plotFr')
axis xy
axis square
colorbar; caxis([0,0.1])
xlabel(xlab);
ylabel(ylab);
title('Standard deviation')

subplot(223);
% plotFr = meanFreqi./repmat(xax',[1,length(yax)]);
plotFr = 1-2*abs(0.5 - meanFreqi.^2 ./(meanFreqi.^2+repmat(xax',[1,length(yax)]).^2));
binFr = zeros(size(plotFr));
binFr(plotFr > 0.95) = 1;
imagesc(xax, yax, plotFr')
axis xy
axis square
hold on
contour(xax, yax' ,binFr',[0.5,0.5], 'w')
colorbar; caxis([0,1])
xlabel(xlab);
ylabel(ylab);
title('Mean')

subplot(224);
% plotFr = stdFreqi./repmat(xax',[1,length(yax)]);
plotFr = 1-2*abs(0.5 - stdFreqi.^2 ./(stdFreqi.^2+repmat(xax',[1,length(yax)]).^2));
imagesc(xax, yax, plotFr')
axis xy
axis square
colorbar; caxis([0,0.1])
xlabel(xlab);
ylabel(ylab);
title('Standard deviation')
