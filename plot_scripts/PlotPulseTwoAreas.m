% pulse plots

%% find phases

Tper = min(floor(1000./avgmorfreq(1, :)));

[~, Id] = findpeaks(phaseOfPulse(2,:,1), 'npeaks', 1);

if Id < Tper+1
    PhaseIds = [Id+1:Tper+1, 1:Id];
else
    PhaseIds = 1:Tper+1;
end

phases = phaseOfPulse(1, PhaseIds');
if phases(1) > phases(2);
    phases = phases(2:end);
    PhaseIds = PhaseIds(2:end)-1;
end % phases(1,1) = NaN; end

% figure;
% plot(phases, '*k ');

%% plots

temp1 = 5:40;
temp2 = temp1(~ismember(temp1, 8:4:40));

% phase differences
minId = 4;
maxId = 40;
phasedata = phasediff(minId:maxId,2);

PhaseIn = wrapTo2Pi((repmat(phases, [length(phasedata),1]) - repmat(phasedata, [1,length(phases)]) + 0.60)*2*pi)/(2*pi);

zax = minId+1:3:maxId-3;
xax = unwrap(repmat(phases,[1,3])*2*pi-2*pi)/(2*pi);
xax1 = unwrap(phases*2*pi)/(2*pi);
colors = parula(length(zax)+2);
% this way, the plot continues to 0 and 1, by adding a period before and
% after; this can be done because the effect is circular. Having the
% figures run from 0 to 1 makes them easier to interpret.



%% significance - region 1

if exist('runnumbers', 'var')
    % diff from 0 - reg 1
    pspike1 = nan(length(PhaseIds),1);
    ptime1 = nan(length(PhaseIds),1);
    pheight1 = nan(length(PhaseIds),1);
    for j = 1:length(PhaseIds)
        pspike1(j) = ttest(h_spikesLRR(1,PhaseIds(j),1,:)) == 1;
        ptime1(j) = ttest(h_peaktimeLRR(1,PhaseIds(j),1,:)) == 1;
        pheight1(j) = ttest(h_peakheightLRR(1,PhaseIds(j),1,:)) == 1;
    end
    
    % spikes
    dat = [];
    adat = [];
    for j = 1:length(PhaseIds)
        dat = [dat, squeeze(h_spikesLRR(1,PhaseIds(j),1,:))'];
        adat = [adat, repmat(xax1(PhaseIds(j)), [1,length(runnumbers)])];
    end
    [rho_spike, p_spike] = circ_corrcl(adat(:), dat(:))
    
    % heights
    dat = [];
    adat = [];
    for j = 1:length(PhaseIds)
        dat = [dat, squeeze(h_peakheightLRR(1,PhaseIds(j),1,:))'];
        adat = [adat, repmat(xax1(PhaseIds(j)), [1,length(runnumbers)])];
    end
    [rho_height, p_height] = circ_corrcl(adat(:), dat(:))
    
    % times
    dat = [];
    adat = [];
    for j = 1:length(PhaseIds)
        dat = [dat, squeeze(h_peaktimeLRR(1,PhaseIds(j),1,:))'];
        adat = [adat, repmat(xax1(PhaseIds(j)), [1,length(runnumbers)])];
    end
    [rho_time, p_time] = circ_corrcl(adat(:), dat(:))
    
end

%% significance - region 2

if exist('runnumbers', 'var')
    % max diff - max - reg 2
    pspikez = nan(length(zax),length(zax));
    pheightzinc = nan(length(zax),length(zax));
    for i = 1:length(zax)
        for j = i+1:length(zax)
            [~,p1] = max(h_spikes(zax(i),PhaseIds,2));
            [~,p2] = max(h_spikes(zax(j),PhaseIds,2));
            
            [~,pspikez(i,j)] = ttest(h_spikesLRR(zax(i),PhaseIds(p1),2,:),h_spikesLRR(zax(j),PhaseIds(p2),2,:));
            [~,pheightzinc(i,j)] = ttest(h_peakheightLRR(zax(i),PhaseIds(p1),2,:),h_peakheightLRR(zax(j),PhaseIds(p2),2,:));
        end
    end
    
    % max diff - min - reg 2
    ptimez = nan(length(zax),length(zax));
    pheightzdec = nan(length(zax),length(zax));
    for i = 1:length(zax)
        for j = i+1:length(zax)
            [~,p1] = min(h_spikes(zax(i),PhaseIds,2));
            [~,p2] = min(h_spikes(zax(j),PhaseIds,2));

            [~,ptimez(i,j)] = ttest(h_peaktimeLRR(zax(i),PhaseIds(p1),2,:),h_peaktimeLRR(zax(j),PhaseIds(p2),2,:));
            [~,pheightzdec(i,j)] = ttest(h_peakheightLRR(zax(i),PhaseIds(p1),2,:),h_peakheightLRR(zax(j),PhaseIds(p2),2,:));
        end
    end
    
    
    % spike - max
    dat = [];
    for i = 1:length(zax)
        [~,pmax] = max(h_spikes(zax(i),PhaseIds,2));
        dat = [dat, squeeze(h_spikesLRR(zax(i),pmax,2,:))];
    end
    
    adat = repmat(phasediff(zax),[length(runnumbers),1]);
    [rho_spike, p_spike] = circ_corrcl(adat(:), dat(:))
    
    % height - max
    dat = [];
    for i = 1:length(zax)
        [~,pmax] = max(h_peakheight(zax(i),PhaseIds,2));
        dat = [dat, squeeze(h_peakheightLRR(zax(i),pmax,2,:))];
    end
    
    adat = repmat(phasediff(zax),[length(runnumbers),1]);
    [rho_heightinc, p_heightinc] = circ_corrcl(adat(:), dat(:))
    
    % height - max
    dat = [];
    for i = 1:length(zax)
        [~,pmin] = min(h_peakheight(zax(i),PhaseIds,2));
        dat = [dat, squeeze(h_peakheightLRR(zax(i),pmin,2,:))];
    end
    
    adat = repmat(phasediff(zax),[length(runnumbers),1]);
    [rho_heightdec, p_heightdec] = circ_corrcl(adat(:), dat(:))
    
    % time - max
    dat = [];
    for i = 1:length(zax)
        [~,pmin] = min(h_peaktime(zax(i),PhaseIds,2));
        dat = [dat, squeeze(h_peaktimeLRR(zax(i),pmin,2,:))];
    end
    
    adat = repmat(phasediff(zax),[length(runnumbers),1]);
    [rho_timedec, p_timedec] = circ_corrcl(adat(:), dat(:))
    
end


%% significance plots

% spikes
[pthr,pcor,padj] = fdr(pspikez(~isnan(pspikez(:))));
dum = find(triu(ones(size(pspikez)),1));
pspikezcor = nan(size(pspikez));
pspikezcor(dum) = padj;

figure;
hold on
de = 0.01;
gray = [0.8,0.8,0.8];
for i = 1:length(zax)
    for j = i+1:length(zax)
        if pspikezcor(i,j) <= 0.05 %/ ((length(zax)*(length(zax)-1))/2)
            plot([phasediff(zax(i),1)-de, phasediff(zax(i),1)+de], [phasediff(zax(j),1),phasediff(zax(j),1)], 'color', colors(j,:), 'linewidth', 3)
            plot([phasediff(zax(i),1),phasediff(zax(i),1)], [phasediff(zax(j),1)-de,phasediff(zax(j),1)+de], 'color', colors(i,:), 'linewidth', 3)
            
            plot([phasediff(zax(j),1)-de,phasediff(zax(j),1)+de],[phasediff(zax(i),1),phasediff(zax(i),1)],  'color', colors(i,:), 'linewidth', 3)
            plot([phasediff(zax(j),1),phasediff(zax(j),1)],[phasediff(zax(i),1)-de, phasediff(zax(i),1)+de],  'color', colors(j,:), 'linewidth', 3)
        else
            plot([phasediff(zax(j),1),phasediff(zax(j),1)],[phasediff(zax(i),1)-de/2, phasediff(zax(i),1)+de/2],  'color', gray, 'linewidth', 1)
            plot([phasediff(zax(j),1)-de/2,phasediff(zax(j),1)+de/2],[phasediff(zax(i),1),phasediff(zax(i),1)],  'color', gray, 'linewidth', 1)
            plot([phasediff(zax(i),1)-de/2, phasediff(zax(i),1)+de/2], [phasediff(zax(j),1),phasediff(zax(j),1)], 'color', gray, 'linewidth', 1)
            plot([phasediff(zax(i),1),phasediff(zax(i),1)], [phasediff(zax(j),1)-de/2,phasediff(zax(j),1)+de/2], 'color', gray, 'linewidth', 1)
        
%             plot(phasediff(zax(i),1), phasediff(zax(j),1), 'color', [0.7,0.7,0.7])
        end
    end
end
plot([0,1],[0,1], 'k:')
xlim([0.3,0.8])
ylim([0.3,0.8])
axis square
xlabel('\Delta\phi')
ylabel('\Delta\phi')
title('Max. \Delta spike rate')


% alternative visualization
% for i = 1:length(zax)
%     [~,p1] = max(h_spikes(zax(i),PhaseIds,2));
%     pp(i) = xax1(p1);
% end
% pp
% pp = [0.6720,0.6077,0.5977,0.5429,0.5329,0.5229,0.4697,0.4597,0.4497,0.3866,0.3145];
% 
% figure;
% hold on
% de = 0.01;
% gray = [0.8,0.8,0.8];
% for i = 1:length(zax)
%     [~,p1] = max(h_spikes(zax(i),PhaseIds,2));
%     plot([pp(i),pp(i)], [0,1.5], 'color', colors(i,:), 'linewidth', 1)
%     plot([0,1], [h_spikes(zax(i),PhaseIds(p1),2),h_spikes(zax(i),PhaseIds(p1),2)], 'color', colors(i,:), 'linewidth', 1)
% end
% for i = 1:length(zax)
%     for j = i+1:length(zax)
%         [~,p1] = max(h_spikes(zax(i),PhaseIds,2));
%         [~,p2] = max(h_spikes(zax(j),PhaseIds,2));
%         if pspikez(i,j) <= 0.05
%             plot(pp(i), h_spikes(zax(j),PhaseIds(p2),2), '*','color', 'k', 'MarkerSize', 8, 'linewidth', 2)
%             plot(pp(j), h_spikes(zax(i),PhaseIds(p1),2), '*','color', 'k','MarkerSize', 8, 'linewidth', 2)
%         
% %             plot([pp(i)-de,pp(i)+de], [h_spikes(zax(j),PhaseIds(p2),2),h_spikes(zax(j),PhaseIds(p2),2)], 'color', colors(j,:), 'linewidth', 3)
% %             plot([pp(j)-de,pp(j)+de], [h_spikes(zax(i),PhaseIds(p1),2),h_spikes(zax(i),PhaseIds(p1),2)], 'color', colors(i,:), 'linewidth', 3)
%         else
%         plot(pp(i), h_spikes(zax(j),PhaseIds(p2),2), '*','color', gray, 'MarkerSize', 8, 'linewidth', 2)
%             plot(pp(j), h_spikes(zax(i),PhaseIds(p1),2), '*','color', gray,'MarkerSize', 8, 'linewidth', 2)    
%         end
%     end
% end
% for i = 1:length(zax)
%     [~,p1] = max(h_spikes(zax(i),PhaseIds,2));
% %     plot([pp(i)-2*de,pp(i)+2*de], [h_spikes(zax(i),PhaseIds(p1),2),h_spikes(zax(i),PhaseIds(p1),2)],'color', colors(i,:), 'linewidth', 3)
%     plot([pp(i)], [h_spikes(zax(i),PhaseIds(p1),2),h_spikes(zax(i),PhaseIds(p1),2)],'s', 'color', colors(i,:), 'linewidth', 4, 'MarkerFaceColor', colors(i,:))
% end

% peakheight - increase
[pthr,pcor,padj] = fdr(pheightzinc(~isnan(pheightzinc(:))));
dum = find(triu(ones(size(pheightzinc)),1));
pheightzinccor = nan(size(pheightzinc));
pheightzinccor(dum) = padj;


figure;
hold on
de = 0.01;
gray = [0.8,0.8,0.8];
for i = 1:length(zax)
    for j = i+1:length(zax)
        if pheightzinccor(i,j) <= 0.05 % / ((length(zax)*(length(zax)-1))/2)
            plot([phasediff(zax(i),1)-de, phasediff(zax(i),1)+de], [phasediff(zax(j),1),phasediff(zax(j),1)], 'color', colors(j,:), 'linewidth', 3)
            plot([phasediff(zax(i),1),phasediff(zax(i),1)], [phasediff(zax(j),1)-de,phasediff(zax(j),1)+de], 'color', colors(i,:), 'linewidth', 3)
            
            plot([phasediff(zax(j),1)-de,phasediff(zax(j),1)+de],[phasediff(zax(i),1),phasediff(zax(i),1)],  'color', colors(i,:), 'linewidth', 3)
            plot([phasediff(zax(j),1),phasediff(zax(j),1)],[phasediff(zax(i),1)-de, phasediff(zax(i),1)+de],  'color', colors(j,:), 'linewidth', 3)
        else
            plot([phasediff(zax(j),1),phasediff(zax(j),1)],[phasediff(zax(i),1)-de/2, phasediff(zax(i),1)+de/2],  'color', gray, 'linewidth', 1)
            plot([phasediff(zax(j),1)-de/2,phasediff(zax(j),1)+de/2],[phasediff(zax(i),1),phasediff(zax(i),1)],  'color', gray, 'linewidth', 1)
            plot([phasediff(zax(i),1)-de/2, phasediff(zax(i),1)+de/2], [phasediff(zax(j),1),phasediff(zax(j),1)], 'color', gray, 'linewidth', 1)
            plot([phasediff(zax(i),1),phasediff(zax(i),1)], [phasediff(zax(j),1)-de/2,phasediff(zax(j),1)+de/2], 'color', gray, 'linewidth', 1)
        
%             plot(phasediff(zax(i),1), phasediff(zax(j),1), 'color', [0.7,0.7,0.7])
        end
    end
end
plot([0,1],[0,1], 'k:')
xlim([0.3,0.8])
ylim([0.3,0.8])
axis square
xlabel('\Delta\phi')
ylabel('\Delta\phi')
title('Max. \Delta peak height')

% peakheight - decrease
[pthr,pcor,padj] = fdr(pheightzdec(~isnan(pheightzdec(:))));
dum = find(triu(ones(size(pheightzdec)),1));
pheightzdeccor = nan(size(pheightzdec));
pheightzdeccor(dum) = padj;

figure;
hold on
de = 0.01;
gray = [0.8,0.8,0.8];
for i = 1:length(zax)
    for j = i+1:length(zax)
        if pheightzdeccor(i,j) <= 0.05
            plot([phasediff(zax(i),1)-de, phasediff(zax(i),1)+de], [phasediff(zax(j),1),phasediff(zax(j),1)], 'color', colors(j,:), 'linewidth', 3)
            plot([phasediff(zax(i),1),phasediff(zax(i),1)], [phasediff(zax(j),1)-de,phasediff(zax(j),1)+de], 'color', colors(i,:), 'linewidth', 3)
            
            plot([phasediff(zax(j),1)-de,phasediff(zax(j),1)+de],[phasediff(zax(i),1),phasediff(zax(i),1)],  'color', colors(i,:), 'linewidth', 3)
            plot([phasediff(zax(j),1),phasediff(zax(j),1)],[phasediff(zax(i),1)-de, phasediff(zax(i),1)+de],  'color', colors(j,:), 'linewidth', 3)
        else
            plot([phasediff(zax(j),1),phasediff(zax(j),1)],[phasediff(zax(i),1)-de/2, phasediff(zax(i),1)+de/2],  'color', gray, 'linewidth', 1)
            plot([phasediff(zax(j),1)-de/2,phasediff(zax(j),1)+de/2],[phasediff(zax(i),1),phasediff(zax(i),1)],  'color', gray, 'linewidth', 1)
            plot([phasediff(zax(i),1)-de/2, phasediff(zax(i),1)+de/2], [phasediff(zax(j),1),phasediff(zax(j),1)], 'color', gray, 'linewidth', 1)
            plot([phasediff(zax(i),1),phasediff(zax(i),1)], [phasediff(zax(j),1)-de/2,phasediff(zax(j),1)+de/2], 'color', gray, 'linewidth', 1)
        
%             plot(phasediff(zax(i),1), phasediff(zax(j),1), 'color', [0.7,0.7,0.7])
        end
    end
end
plot([0,1],[0,1], 'k:')
xlim([0.3,0.8])
ylim([0.3,0.8])
axis square
xlabel('\Delta\phi')
ylabel('\Delta\phi')
title('Min. \Delta peak height')


[pthr,pcor,padj] = fdr(ptimez(~isnan(ptimez(:))));
dum = find(triu(ones(size(ptimez)),1));
ptimezcor = nan(size(ptimez));
ptimezcor(dum) = padj;

figure;
hold on
de = 0.01;
gray = [0.8,0.8,0.8];
for i = 1:length(zax)
    for j = i+1:length(zax)
        if ptimezcor(i,j) <= 0.05
            plot([phasediff(zax(i),1)-de, phasediff(zax(i),1)+de], [phasediff(zax(j),1),phasediff(zax(j),1)], 'color', colors(j,:), 'linewidth', 3)
            plot([phasediff(zax(i),1),phasediff(zax(i),1)], [phasediff(zax(j),1)-de,phasediff(zax(j),1)+de], 'color', colors(i,:), 'linewidth', 3)
            
            plot([phasediff(zax(j),1)-de,phasediff(zax(j),1)+de],[phasediff(zax(i),1),phasediff(zax(i),1)],  'color', colors(i,:), 'linewidth', 3)
            plot([phasediff(zax(j),1),phasediff(zax(j),1)],[phasediff(zax(i),1)-de, phasediff(zax(i),1)+de],  'color', colors(j,:), 'linewidth', 3)
        else
            plot([phasediff(zax(j),1),phasediff(zax(j),1)],[phasediff(zax(i),1)-de/2, phasediff(zax(i),1)+de/2],  'color', gray, 'linewidth', 1)
            plot([phasediff(zax(j),1)-de/2,phasediff(zax(j),1)+de/2],[phasediff(zax(i),1),phasediff(zax(i),1)],  'color', gray, 'linewidth', 1)
            plot([phasediff(zax(i),1)-de/2, phasediff(zax(i),1)+de/2], [phasediff(zax(j),1),phasediff(zax(j),1)], 'color', gray, 'linewidth', 1)
            plot([phasediff(zax(i),1),phasediff(zax(i),1)], [phasediff(zax(j),1)-de/2,phasediff(zax(j),1)+de/2], 'color', gray, 'linewidth', 1)
        
%             plot(phasediff(zax(i),1), phasediff(zax(j),1), 'color', [0.7,0.7,0.7])
        end
    end
end
plot([0,1],[0,1], 'k:')
xlim([0.3,0.8])
ylim([0.3,0.8])
axis square
xlabel('\Delta\phi')
ylabel('\Delta\phi')
title('Min. \Delta peak times')


%% line plots

% h_spikes
figure;
subplot(311)
plot(xax, repmat(h_spikes(1,PhaseIds,1), [1,3,1]), 'color', [0.,0.,0.], 'linewidth', 2);
hold on
errorbar(xax, repmat(h_spikes(1,PhaseIds,1), [1,3,1]), repmat(h_spikes_std(1,PhaseIds,1), [1,3,1]), 'color', [0.,0.,0.], 'linewidth', 1);
psign = xax;
psign(repmat(pspike1,[1,3])==0) = NaN;
plot(psign, 0*ones(length(xax),1),'k', 'linewidth',3);
ylim([0,1])
xlim([0,1])
xlabel('\theta_{1}')

subplot(312);
hold on;
for i = 1:length(zax)
    plot(xax, repmat(h_spikes(zax(i),PhaseIds,2), [1,3,1]), 'color', colors(i,:), 'linewidth', 2);
end
ylim([0,1.5])
xlim([0,1])
ylabel('Change in spike rate (Hz)')
xlabel('\theta_{1}')
legend(strcat('\Delta\phi = ', num2str(phasediff(zax,1), 2)))

subplot(313);
hold on;
for i = 1:length(zax)
    [phase2, IdsSort] = sort(wrapTo2Pi((phases - phasediff(zax(i),1))*2*pi) / (2*pi));
    plot(unwrap(repmat(phase2,[1,3])*2*pi-2*pi)/(2*pi), repmat(h_spikes(zax(i),PhaseIds(IdsSort),2), [1,3,1]), 'color', colors(i,:), 'linewidth', 2);
end
ylim([0,1.5])
xlim([0,1])
xlabel('\theta_{2}')

% h_peak
% figure;
% subplot(311)
% plot(xax, repmat(h_peak(1,PhaseIds,1), [1,3,1]), 'color', [0.,0.,0.], 'linewidth', 2);
% hold on
% errorbar(xax, repmat(h_peak(1,PhaseIds,1), [1,3,1]), repmat(h_peak_std(1,PhaseIds,1), [1,3,1]), 'color', [0.,0.,0.], 'linewidth', 1);
% ylim([-50,50])
% xlim([0,1])
% xlabel('\theta_{1}')
%
% subplot(312);
% hold on;
% for i = 1:length(zax)
%     plot(xax, repmat(h_peak(zax(i),PhaseIds,2), [1,3,1]), 'color', colors(i,:), 'linewidth', 2);
% end
% ylim([-50,50])
% xlim([0,1])
% ylabel('Change in peak height (ms)')
% xlabel('\theta_{1}')
% legend(strcat('\Delta\phi = ', num2str(phasediff(zax,1), 2)))
%
% subplot(313);
% hold on;
% for i = 1:length(zax)
%     [phase2, IdsSort] = sort(wrapTo2Pi((phases - phasediff(zax(i),1))*2*pi) / (2*pi));
%     plot(unwrap(repmat(phase2,[1,3])*2*pi-2*pi)/(2*pi), repmat(h_peak(zax(i),PhaseIds(IdsSort),2), [1,3,1]), 'color', colors(i,:), 'linewidth', 2);
% end
% ylim([-50,50])
% xlim([0,1])
% xlabel('\theta_{2}')

% h_peakheight
figure;
subplot(311)
plot(xax, repmat(h_peakheight(1,PhaseIds,1), [1,3,1]), 'color', [0.,0.,0.], 'linewidth', 2);
hold on
errorbar(xax, repmat(h_peakheight(1,PhaseIds,1), [1,3,1]), repmat(h_peakheight_std(1,PhaseIds,1), [1,3,1]), 'color', [0.,0.,0.], 'linewidth', 1);
psign = xax;
psign(repmat(pheight1,[1,3])==0) = NaN;
plot(psign, -100*ones(length(xax),1),'k', 'linewidth',3);
ylim([-100,100])
xlim([0,1])
xlabel('\theta_{1}')

subplot(312);
hold on;
for i = 1:length(zax)
    plot(xax, repmat(h_peakheight(zax(i),PhaseIds,2), [1,3,1]), 'color', colors(i,:), 'linewidth', 2);
end
ylim([-100,100])
xlim([0,1])
ylabel('Change in peak height (spikes/s)')
xlabel('\theta_{1}')
legend(strcat('\Delta\phi = ', num2str(phasediff(zax,1), 2)))

subplot(313);
hold on;
for i = 1:length(zax)
    [phase2, IdsSort] = sort(wrapTo2Pi((phases - phasediff(zax(i),1))*2*pi) / (2*pi));
    plot(unwrap(repmat(phase2,[1,3])*2*pi-2*pi)/(2*pi), repmat(h_peakheight(zax(i),PhaseIds(IdsSort),2), [1,3,1]), 'color', colors(i,:), 'linewidth', 2);
end
ylim([-100,100])
xlim([0,1])
xlabel('\theta_{2}')

% h_power
% figure;
% subplot(311)
% plot(xax, repmat(h_power(1,PhaseIds,1), [1,3,1]), 'color', [0.6,0.6,0.6], 'linewidth', 2);
% hold on
% errorbar(xax, repmat(h_power(1,PhaseIds,1), [1,3,1]), repmat(h_power_std(1,PhaseIds,1), [1,3,1]), 'color', [0.6,0.6,0.6], 'linewidth', 2);
% ylim([-100,50])
% xlim([0,1])
% xlabel('\theta_{1}')
%
% subplot(312);
% hold on;
% for i = 1:length(zax)
%     plot(xax, repmat(h_power(zax(i),PhaseIds,2), [1,3,1]), 'color', colors(i,:), 'linewidth', 2);
% end
% ylim([-100,50])
% xlim([0,1])
% ylabel('Change in power')
% xlabel('\theta_{1}')
% legend(strcat('\Delta\phi = ', num2str(phasediff(zax,1), 2)))
%
% subplot(313);
% hold on;
% for i = 1:length(zax)
%     [phase2, IdsSort] = sort(wrapTo2Pi((phases - phasediff(zax(i),1))*2*pi) / (2*pi));
%     plot(unwrap(repmat(phase2,[1,3])*2*pi-2*pi)/(2*pi), repmat(h_power(zax(i),PhaseIds(IdsSort),2), [1,3,1]), 'color', colors(i,:), 'linewidth', 2);
% end
% ylim([-100,50])
% xlim([0,1])
% xlabel('\theta_{2}')


% h_peaktime
figure;
subplot(311)
plot(xax, repmat(h_peaktime(1,PhaseIds,1), [1,3,1]), 'color', [0.,0.,0.], 'linewidth', 2);
hold on
errorbar(xax, repmat(h_peaktime(1,PhaseIds,1), [1,3,1]), repmat(h_peaktime_std(1,PhaseIds,1), [1,3,1]), 'color', [0.,0.,0.], 'linewidth', 1);
psign = xax;
psign(repmat(ptime1,[1,3])==0) = NaN;
plot(psign, -3*ones(length(xax),1),'k', 'linewidth',3);
ylim([-3,1])
xlim([0,1])
xlabel('\theta_{1}')

subplot(312);
hold on;
for i = 1:length(zax)
    plot(xax, repmat(h_peaktime(zax(i),PhaseIds,2), [1,3,1]), 'color', colors(i,:), 'linewidth', 2);
end
ylim([-3,1])
xlim([0,1])
ylabel('Change in peak time (ms)')
xlabel('\theta_{1}')
legend(strcat('\Delta\phi = ', num2str(phasediff(zax,1), 2)))

subplot(313);
hold on;
for i = 1:length(zax)
    [phase2, IdsSort] = sort(wrapTo2Pi((phases - phasediff(zax(i),1))*2*pi) / (2*pi));
    plot(unwrap(repmat(phase2,[1,3])*2*pi-2*pi)/(2*pi), repmat(h_peaktime(zax(i),PhaseIds(IdsSort),2), [1,3,1]), 'color', colors(i,:), 'linewidth', 2);
end
ylim([-3,1])
xlim([0,1])
xlabel('\theta_{2}')


% h_freq
% figure;
% subplot(311)
% plot(xax, repmat(h_freq(1,PhaseIds,1), [1,3,1]), 'color', [0.6,0.6,0.6], 'linewidth', 2);
% hold on
% errorbar(xax, repmat(h_freq(1,PhaseIds,1), [1,3,1]), repmat(h_freq_std(1,PhaseIds,1), [1,3,1]), 'color', [0.6,0.6,0.6], 'linewidth', 2);
% ylim([0,5])
% xlim([0,1])
% xlabel('\theta_{1}')
%
% subplot(312);
% hold on;
% for i = 1:length(zax)
%     plot(xax, repmat(h_freq(zax(i),PhaseIds,2), [1,3,1]), 'color', colors(i,:), 'linewidth', 2);
% end
% ylim([0,5])
% xlim([0,1])
% ylabel('Change in frequecy (Hz)')
% xlabel('\theta_{1}')
% legend(strcat('\Delta\phi = ', num2str(phasediff(zax,1), 2)))
%
% subplot(313);
% hold on;
% for i = 1:length(zax)
%     [phase2, IdsSort] = sort(wrapTo2Pi((phases - phasediff(zax(i),1))*2*pi) / (2*pi));
%     plot(unwrap(repmat(phase2,[1,3])*2*pi-2*pi)/(2*pi), repmat(h_freq(zax(i),PhaseIds(IdsSort),2), [1,3,1]), 'color', colors(i,:), 'linewidth', 2);
% end
% ylim([0,5])
% xlim([0,1])
% xlabel('\theta_{2}')

% h_sync
% figure;
% subplot(311)
% plot(xax, repmat(h_sync(1,PhaseIds,1), [1,3,1]), 'color', [0.6,0.6,0.6], 'linewidth', 2);
% hold on
% errorbar(xax, repmat(h_sync(1,PhaseIds,1), [1,3,1]), repmat(h_sync_std(1,PhaseIds,1), [1,3,1]), 'color', [0.6,0.6,0.6], 'linewidth', 2);
% ylim([-0.3,0.3])
% xlim([0,1])
% xlabel('\theta_{1}')
%
% subplot(312);
% hold on;
% for i = 1:length(zax)
%     plot(xax, repmat(h_sync(zax(i),PhaseIds,2), [1,3,1]), 'color', colors(i,:), 'linewidth', 2);
% end
% ylim([-0.3,0.3])
% xlim([0,1])
% ylabel('Change in synchrony (Hz)')
% xlabel('\theta_{1}')
% legend(strcat('\Delta\phi = ', num2str(phasediff(zax,1), 2)))
%
% subplot(313);
% hold on;
% for i = 1:length(zax)
%     [phase2, IdsSort] = sort(wrapTo2Pi((phases - phasediff(zax(i),1))*2*pi) / (2*pi));
%     plot(unwrap(repmat(phase2,[1,3])*2*pi-2*pi)/(2*pi), repmat(h_sync(zax(i),PhaseIds(IdsSort),2), [1,3,1]), 'color', colors(i,:), 'linewidth', 2);
% end
% ylim([-0.3,0.3])
% xlim([0,1])
% xlabel('\theta_{2}')

%% imagesc


% figure('Position', [50,50,400,600]);
% subplot(10,1,1); imagesc(phases, 1, phaseOfPulse(1,PhaseIds)); set(gca,'ydir', 'normal'); colorbar; caxis([0,1])
% colormap(parula)
% title('Phase of pulse relative to area 2')
% subplot(10,1,2:10); imagesc(phases, phasedata, PhaseIn(:, PhaseIds, 1)); set(gca,'ydir', 'normal'); colorbar; caxis([0,1])
% colormap(parula)
% xlabel('Phase of pulse'); ylabel('Phase difference')
%
% figure('Position', [50,50,400,600]);
% subplot(10,1,1); imagesc(phases, 1, h_spikes(1,PhaseIds,1)); set(gca,'ydir', 'normal'); colorbar; caxis([0,500])
% title('Spikes in first period after pulse');
% subplot(10,1,2:10); imagesc(phases, phasediff(:,1), h_spikes(:,PhaseIds,2)); set(gca,'ydir', 'normal'); colorbar; caxis([0,500])
% xlabel('Phase of pulse'); ylabel('Phase difference')
%
% figure('Position', [50,50,400,600]);
% subplot(10,4,1:3); imagesc(phases, 1, h_spikes(1,PhaseIds,1)); set(gca,'ydir', 'normal');
% colorbar; caxis([0,200])
% colormap(parula)
% title('Period spikes');
% subplot(10,4,temp2); pcolor(phases, phasedata, h_spikes(minId:maxId,PhaseIds,2)); set(gca,'ydir', 'normal');
% colorbar; caxis([0,500])
% xlabel('Phase of pulse'); ylabel('Phase difference')
% shading flat
% colormap(parula)
% subplot(10,4,8:4:40); plot(sum(h_spikes(minId:maxId,PhaseIds,2),2), phasedata, 'color', sett.purple, 'linewidth', 2);
% colormap(parula)
% %xlim([0.5,1.75]); ylim([phasedata(end,1), phasedata(1,1)]);
%
%
% figure('Position', [50,50,400,600]);
% subplot(10,4,1:3); imagesc(phases, 1, h_base(1,PhaseIds,1)); set(gca,'ydir', 'normal'); colorbar; caxis([0,8])
% colormap(parula)
% title('Baseline firing');
% subplot(10,4,temp2); pcolor(phases, phasedata, h_base(minId:maxId,PhaseIds,2)); set(gca,'ydir', 'normal'); colorbar; caxis([0,3])
% xlabel('Phase of pulse'); ylabel('Phase difference')
% shading flat
% colormap(parula)
% subplot(10,4,8:4:40); plot(mean(h_base(minId:maxId,PhaseIds,2),2), phasedata, 'color', sett.purple, 'linewidth', 2);
% colormap(parula)
% xlim([0.5,1.75]); ylim([phasedata(end,1), phasedata(1,1)]);
%
% % figure('Position', [50,50,400,600]);
% % subplot(10,1,1); imagesc(phases, 1, h_fper(1,PhaseIds,1));set(gca,'ydir', 'normal'); colorbar; caxis([0,600])
% % title('Highest concentration of spikes');
% % subplot(10,1,2:10); imagesc(phases, phasediff(:,1), h_fper(:,PhaseIds,2)); set(gca,'ydir', 'normal'); colorbar; caxis([0,600])
% % xlabel('Phase of pulse'); ylabel('Phase difference')
%
% figure('Position', [50,50,400,600]);
% subplot(10,4,1:3); imagesc(phases, 1, h_freq(1,PhaseIds,1));set(gca,'ydir', 'normal'); colorbar; caxis([0,4])
% colormap(parula)
% title('Frequency change');
% subplot(10,4,temp2); pcolor(phases, phasedata, h_freq(minId:maxId,PhaseIds,2)); set(gca,'ydir', 'normal'); colorbar; caxis([0,4])
% xlabel('Phase of pulse'); ylabel('Phase difference')
% shading flat
% colormap(parula)
% subplot(10,4,8:4:40); plot(mean(h_freq(minId:maxId,PhaseIds,2),2), phasedata, 'color', sett.blue, 'linewidth', 2);
% xlim([0,5]); ylim([phasedata(end,1), phasedata(1,1)]);
% colormap(parula)
%
% % figure('Position', [50,50,400,600]);
% % subplot(10,1,1); imagesc(phases, 1, h_power(1,PhaseIds,1));set(gca,'ydir', 'normal'); colorbar; caxis([-50,30])
% % title('Power change');
% % subplot(10,1,2:10); imagesc(phases, phasediff(:,1), h_power(:,PhaseIds,2)); set(gca,'ydir', 'normal'); colorbar; caxis([-50,30])
% % xlabel('Phase of pulse'); ylabel('Phase difference')
%
% figure('Position', [50,50,400,600]);
% subplot(10,4,1:3); imagesc(phases, 1, h_peakheight(1,PhaseIds,1));set(gca,'ydir', 'normal'); colorbar; caxis([-50,50])
% colormap(parula)
% title('LFP amplitude');
% subplot(10,4,temp2); pcolor(phases, phasedata, h_peakheight(minId:maxId,PhaseIds,2)); set(gca,'ydir', 'normal'); colorbar; caxis([-50,50])
% xlabel('Phase of pulse'); ylabel('Phase difference')
% shading flat
% colormap(parula)
% subplot(10,4,8:4:40); plot(mean(h_peakheight(minId:maxId,PhaseIds,2),2), phasedata, 'color', sett.green, 'linewidth', 2);
% xlim([-50,50]); ylim([phasedata(end,1), phasedata(1,1)]);
% colormap(parula)
%
% figure('Position', [50,50,400,600]);
% subplot(10,4,1:3); imagesc(phases, 1, h_peaktime(1,PhaseIds,1));
% set(gca,'ydir', 'normal'); colorbar; caxis([-13,13])
% colormap(parula)
% title('LFP amplitude');
% subplot(10,4,temp2); pcolor(phases, phasedata, h_peaktime(minId:maxId,PhaseIds,2));
% set(gca,'ydir', 'normal'); colorbar; caxis([-13,13])
% xlabel('Phase of pulse'); ylabel('Phase difference')
% shading flat
% colormap(parula)
% subplot(10,4,8:4:40); plot(mean(h_peaktime(minId:maxId,PhaseIds,2),2), phasedata, 'color', sett.green, 'linewidth', 2);
% xlim([-13,13]); ylim([phasedata(end,1), phasedata(1,1)]);
% colormap(parula)
%
% figure('Position', [50,50,400,600]);
% subplot(10,4,1:3); imagesc(phases, 1, h_sync(1,PhaseIds,1));set(gca,'ydir', 'normal');
% colorbar; caxis([-0.5,0.5])
% colormap(parula)
% title('LFP amplitude');
% subplot(10,4,temp2); pcolor(phases, phasedata, h_sync(minId:maxId,PhaseIds,2));
% set(gca,'ydir', 'normal'); colorbar; caxis([-0.5,0.5])
% xlabel('Phase of pulse'); ylabel('Phase difference')
% shading flat
% colormap(parula)
% subplot(10,4,8:4:40); plot(mean(h_sync(minId:maxId,PhaseIds,2),2), phasedata, 'color', sett.green, 'linewidth', 2);
% xlim([-0.5,0.5]); ylim([phasedata(end,1), phasedata(1,1)]);
% colormap(parula)
%
% % figure('Position', [50,50,400,600]);
% % subplot(10,1,1); imagesc(phases, 1, CVe(1,PhaseIds,1)); set(gca,'ydir', 'normal'); colorbar; caxis([0.04,0.16])
% % title('CV')
% % subplot(10,1,2:10); imagesc(phases, phasediff(:,1), CVe(:,PhaseIds,2)); set(gca,'ydir', 'normal'); colorbar; caxis([0.04,0.16])
% % xlabel('Phase of pulse'); ylabel('Phase difference')
%
