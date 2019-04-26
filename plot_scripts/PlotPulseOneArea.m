% pulse plots

%% find phases

if sett.Iapploop >= 2
    Ids = 1:length(Ids);
    ids = Ids; %[Ids,(Ids(end)+1)];
    phases = phaseOfPulse(ids,1)';
end
% if phases(1,1) > 0.5; phases(1,1) = NaN; end

% figure;
% plot(phases, '*k ');

%% for different pulse amp

if sett.Iapploop >= 2
    
    zax = [5:4:41];
    xax = unwrap(repmat(phases,[1,3])*2*pi-2*pi)/pi;
    colors = flipud(gray(length(zax)+5));
    % this way, the plot continues to 0 and 1, by adding a period before and
    % after; this can be done because the effect is circular. Having the
    % figures run from 0 to 1 makes them easier to interpret.
    
    % h_spikes
    figure;
    hold on;
    for i = 1:length(zax)
        plot(xax, repmat(h_spikes(ids,zax(i))/5*400, [3,1]), 'color', colors(i+2,:), 'linewidth', 2);
    end
    ylim([-20,120]) %ylim([-0.25,1.5])
    xlim([0,2])
    ylabel('Change in spike rate (Hz)')
    xlabel('\theta_{1}')
    
    
    % h_peakheight
    figure;
    hold on;
    for i = 1:length(zax)
        plot(xax, repmat(h_peakheight(ids,zax(i)), [3,1]), 'color', colors(i+2,:), 'linewidth', 2);
    end
    ylim([-50,150])
    %     xlim([0,2])
    ylabel('Change in peak height (spikes/s)')
    xlabel('\theta_{1}')
    legend(strcat('% stimulated = ', num2str(zax'*2.5-2.5)), 'Location','NorthWest')
    
    
    % h_peaktime
    figure;
    hold on;
    for i = 1:length(zax)
        plot(xax, repmat(h_peaktime(ids,zax(i)), [3,1]), 'color', colors(i+2,:), 'linewidth', 2);
    end
    ylim([-5,1])
    xlim([0,2])
    ylabel('Change in peak time (ms)')
    xlabel('\theta_{1}')
    
end


%% for different synchronizations

if sett.Iapploop == 1 && ~isempty(find(sett.loopIe == 5))
    % rows = sync; columns = phase of pulse;
    
    % zax = [3,6,9,11,12,14,17,21,30]; % rows
    % zax = [1,3,4,5,6,7,10,11,12,13,14,15,16,17,20,21,24,30];
    if reg == 1
        zax = [2,4,5,10,11,12,14,16,20,30];
        colors = flipud(gray(length(zax)+4));
    else
        zax = 1;
        colors = zeros(3,3);
    end
    
    avgppc = mean(Phaselocke(:,:,1,:),2);
    avgre = mean(re(:,:,1,:),2);
    avgfreq = mean(avgmorfreq(:,:,1,:),2);
    
    
    %     figure(1); hold on;
    figure(2); hold on;
    figure(3); hold on;
    figure(4); hold on;
    
    for i = 1:length(zax)
        % find phases
        %                 Tper = floor(1000./avgmorfreq(zax(i), 1));
        [~, Id] = findpeaks(phaseOfPulse(zax(i),:), 'npeaks', 1);
        
        if Id < Tper+1
            Ids = [Id+1:Tper+1, 1:Id];
        else
            Ids = 1:Tper+1;
        end
        
        ids = Ids;
        phases = phaseOfPulse(zax(i),ids)';
        if phases(1,1) >= phases(2,1);%0.5;
            phases = phases(2:end);
            ids = Ids(2:end);
        end
        
        if exist('runnumbers', 'var')
            % stats - reg 1
            for j = 1:length(ids)
                [~,pspike(j)] = ttest(h_spikesLRR(zax(i),ids(j),1,:));
                [~,ptime(j)] = ttest(h_peaktimeLRR(zax(i),ids(j),1,:));
                [~,pheight(j)] = ttest(h_peakheightLRR(zax(i),ids(j),1,:));
%                 ids(i)
%                 [~,pspike(j)] = ttest(squeeze(h_spikes0LR(zax(1),ids(j),1,:)), squeeze(h_spikes1LR(zax(1),ids(j),1,:)));
%                 [~,pheight(j)] = ttest(squeeze(h_peakheight0LR(zax(1),ids(j),1,:)), squeeze(h_peakheight1LR(zax(1),ids(j),1,:)));
%                 [~,ptime(j)] = ttest(squeeze(h_peaktime0LR(zax(1),ids(j),1,:)), squeeze(h_peaktime1LR(zax(1),ids(j),1,:)));
            end
        end
        
        xax = unwrap(repmat(phases,[3,1])*2*pi-2*pi)/pi;
        
        %         figure(1); plot(xax, repmat(h_spikes(zax(i), ids)*Ne*dthist, [1,3]), 'color', colors(i+2,:), 'linewidth', 2);
        
        figure(4);
        plot(xax, repmat((h_spikes(zax(i), ids)/(avgre(zax(i))/max(avgre)))/5 *400, [1,3]), 'color', colors(i+2,:), 'linewidth', 2);
        if reg >1
            errorbar(xax, repmat((h_spikes(zax(i), ids)/(avgre(zax(i))/max(avgre)))/5 *400, [1,3]), ...
                repmat((h_spikes_std(zax(i), ids)/(avgre(zax(i))/max(avgre)))/5 *400, [1,3]), 'color', colors(i+2,:), 'linewidth', 1);
            plot(unwrap(repmat(phases(pspike<0.05),[3,1])*2*pi-2*pi)/pi,...
                repmat((h_spikes(zax(i), ids(pspike<0.05))/(avgre(zax(i))/max(avgre)))/5 *400, [1,3]) + ...
                repmat((h_spikes_std(zax(i), ids(pspike<0.05))/(avgre(zax(i))/max(avgre)))/5 *400, [1,3]) + 5,...
                '*', 'color', colors(i+2,:))
        end
        
        figure(2); plot(xax, repmat(h_peakheight(zax(i), ids), [1,3]), 'color', colors(i+2,:), 'linewidth', 2);
        if reg >1
            errorbar(xax, repmat(h_peakheight(zax(i), ids), [1,3]), ...
                repmat(h_peakheight_std(zax(i), ids), [1,3]), 'color', colors(i+2,:), 'linewidth', 1);
            plot(unwrap(repmat(phases(pheight<0.05),[3,1])*2*pi-2*pi)/pi,...
                repmat(h_peakheight(zax(i), ids(pheight<0.05)), [1,3]) + ...
                repmat(h_peakheight_std(zax(i), ids(pheight<0.05)), [1,3]) + 5,...
                '*', 'color', colors(i+2,:))
        end
        
        figure(3); plot(xax, repmat(h_peaktime(zax(i), ids), [1,3]), 'color', colors(i+2,:), 'linewidth', 2);
        if reg >1
            errorbar(xax, repmat(h_peaktime(zax(i), ids), [1,3]), ...
                repmat(h_peaktime_std(zax(i), ids), [1,3]), 'color', colors(i+2,:), 'linewidth', 1);
            plot(unwrap(repmat(phases(ptime<0.05),[3,1])*2*pi-2*pi)/pi,...
                repmat(h_peaktime(zax(i), ids(ptime<0.05)), [1,3]) + ...
                repmat(h_peaktime_std(zax(i), ids(ptime<0.05)), [1,3]) + 0.15,...
                '*', 'color', colors(i+2,:))
        end
        
    end
    
    %     figure(1)
    %     ylim([-20,120])
    %     xlim([0,2])
    %     ylabel('Change in spikes/T')
    %     xlabel('\theta_{1}')
    
    figure(2)
    ylim([-50,100]); %ylim([-50,150])
    xlim([0,2])
    ylabel('Change in peak height (spikes/s)')
    xlabel('\theta_{1}')
    legend(strcat('PPC = ', num2str(avgppc(zax))), 'Location','NorthWest')
    
    figure(3)
    ylim([-3,1]);% ylim([-8,8])
    xlim([0,2])
    ylabel('Change in peak time (ms)')
    xlabel('\theta_{1}')
    
    figure(4)
    ylim([-20,120]);% ylim([-20,120])
    xlim([0,2])
    ylabel('Change in spikes/T (corrected for avg. spike rate)')
    xlabel('\theta_{1}')
    
end

%% for osc. inhibition

if sett.Iapploop == 1 && sett.Iperloop>0
    
    % rows = phase of pulse; columns = amp or per of oscillation
    
    % zax = [3,6,9,11,12,14,17,21,30]; % rows
    % zax = [1,3,4,5,6,7,10,11,12,13,14,15,16,17,20,21,24,30];
    if reg == 1
        zax = 1:2:l2steps;% [2,4,5,10,11,12,14,16,20,30];
        colors = flipud(gray(length(zax)+4));
    else
        zax = 1;
        colors = zeros(3,3);
    end
    
    avgppc = mean(Phaselocke(:,:,1,:),1);
    avgre = mean(re(:,:,1,:),1);
    avgfreq = mean(avgmorfreq(:,:,1,:),1);
    
    
    %     figure(1); hold on;
    figure(2); hold on;
    figure(3); hold on;
    figure(4); hold on;
    
    for i = 1:length(zax)
        % find phases
        %                 Tper = floor(1000./avgmorfreq(zax(i), 1));
        [~, Id] = findpeaks(phaseOfPulse(:,zax(i)), 'npeaks', 1);
        
        if Id < Tper+1
            Ids = [Id+1:Tper+1, 1:Id];
        else
            Ids = 1:Tper+1;
        end
        
        ids = Ids;
        phases = phaseOfPulse(ids,zax(i))';
        if phases(1) >= phases(2);%0.5;
            phases = phases(2:end);
            ids = Ids(2:end);
        end
        
        if exist('runnumbers', 'var')
            % stats - reg 1
            for j = 1:length(ids)
                [~,pspike(j)] = ttest(h_spikesLRR(ids(j),zax(i),1,:));
                [~,ptime(j)] = ttest(h_peaktimeLRR(ids(j),zax(i),1,:));
                [~,pheight(j)] = ttest(h_peakheightLRR(ids(j),zax(i),1,:));
%                 ids(i)
%                 [~,pspike(j)] = ttest(squeeze(h_spikes0LR(zax(1),ids(j),1,:)), squeeze(h_spikes1LR(zax(1),ids(j),1,:)));
%                 [~,pheight(j)] = ttest(squeeze(h_peakheight0LR(zax(1),ids(j),1,:)), squeeze(h_peakheight1LR(zax(1),ids(j),1,:)));
%                 [~,ptime(j)] = ttest(squeeze(h_peaktime0LR(zax(1),ids(j),1,:)), squeeze(h_peaktime1LR(zax(1),ids(j),1,:)));
            end
        end
        
        xax = unwrap(repmat(phases,[1,3])*2*pi-2*pi)/pi;
        
        %         figure(1); plot(xax, repmat(h_spikes(zax(i), ids)*Ne*dthist, [1,3]), 'color', colors(i+2,:), 'linewidth', 2);
        
        figure(4);
        histfact = (Ne*dthist)/1000; % to compensate for the scaling done in the histogram function
        plot(xax, repmat((h_spikes(ids,zax(i))/(avgre(zax(i))/max(avgre)))*histfact, [3,1]), '-o', 'color', colors(i+2,:), 'linewidth', 2);
        if reg >1
            errorbar(xax, repmat((h_spikes(ids, zax(i))/(avgre(zax(i))/max(avgre)))*histfact, [3,1]), ...
                repmat((h_spikes_std(ids, zax(i))/(avgre(zax(i))/max(avgre)))*histfact, [3,1]), 'color', colors(i+2,:), 'linewidth', 1);
            plot(unwrap(repmat(phases(pspike<0.05),[1,3])*2*pi-2*pi)/pi,...
                repmat((h_spikes(ids(pspike<0.05),zax(i))/(avgre(zax(i))/max(avgre)))*histfact , [3,1]) + ...
                repmat((h_spikes_std(ids(pspike<0.05),zax(i))/(avgre(zax(i))/max(avgre)))*histfact, [3,1]) + 5,...
                '*', 'color', colors(i+2,:))
        end
        
        figure(2); plot(xax, repmat(h_peakheight(ids, zax(i)), [3,1]), 'color', colors(i+2,:), 'linewidth', 2);
        if reg >1
            errorbar(xax, repmat(h_peakheight(ids, zax(i)), [3,1]), ...
                repmat(h_peakheight_std(ids, zax(i)), [3,1]), 'color', colors(i+2,:), 'linewidth', 1);
            plot(unwrap(repmat(phases(pheight<0.05),[1,3])*2*pi-2*pi)/pi,...
                repmat(h_peakheight(ids(pheight<0.05),zax(i)), [3,1]) + ...
                repmat(h_peakheight_std(ids(pheight<0.05),zax(i)), [3,1]) + 5,...
                '*', 'color', colors(i+2,:))
        end
        
        figure(3); plot(xax, repmat(h_peaktime(ids, zax(i)), [3,1]), 'color', colors(i+2,:), 'linewidth', 2);
        if reg >1
            errorbar(xax, repmat(h_peaktime(ids, zax(i)), [3,1]), ...
                repmat(h_peaktime_std(ids, zax(i)), [3,1]), 'color', colors(i+2,:), 'linewidth', 1);
            plot(unwrap(repmat(phases(ptime<0.05),[1,3])*2*pi-2*pi)/pi,...
                repmat(h_peaktime(ids(ptime<0.05),zax(i)), [3,1]) + ...
                repmat(h_peaktime_std(ids(ptime<0.05),zax(i)), [3,1]) + 0.15,...
                '*', 'color', colors(i+2,:))
        end
        
    end
    
    %     figure(1)
    %     ylim([-20,120])
    %     xlim([0,2])
    %     ylabel('Change in spikes/T')
    %     xlabel('\theta_{1}')
    
    figure(2)
%     ylim([-50,100]); %ylim([-50,150])
    xlim([0,2])
    ylabel('Change in peak height (spikes/s)')
    xlabel('\theta_{1}')
    legend(strcat('PPC = ', num2str(avgppc(zax))), 'Location','NorthWest')
    
    figure(3)
%     ylim([-3,1]);% ylim([-8,8])
    xlim([0,2])
    ylabel('Change in peak time (ms)')
    xlabel('\theta_{1}')
    
    figure(4)
%     ylim([-20,120]);% ylim([-20,120])
    xlim([0,2])
    ylabel('Change in spikes/T (corrected for avg. spike rate)')
    xlabel('\theta_{1}')
%     
end

%% example traces - for LoadRuns
% xl = [300,450];
% exrun = 1;
%
% for ps = 1:size(phaseOfPulse,2);
%     pulse_start = sett.Ipulse_start + (ps-1)*sett.timestep;
%
%     figure
%
%     subplot(411); hold on
%     plot(dthist:dthist:sett.Ttot, squeeze(histe0_saveLR(1,1,:,exrun)), 'color', [0.6, 0.6,0.6], 'linewidth', 2)
%     plot(dthist:dthist:sett.Ttot, squeeze(histe_saveLR(1,ps,:,exrun)), 'color', sett.red, 'linewidth', 2)
%     plot([pulse_start, pulse_start], [0,400], 'color', [243,161,0]/255, 'linewidth', 2)
%     ylabel('Spike density (spikes/s)')
%     xlim(xl)
%     title(['\theta = ', num2str(phaseOfPulseLR(1,ps,1,exrun))])
%
%     subplot(412); hold on
%     plot(Tselection:dthist:sett.Ttot-50, squeeze(sigfper0_saveLR(1,1,:,exrun))/Ne, 'color', [0.6, 0.6,0.6], 'linewidth', 2)
%     plot(Tselection:dthist:sett.Ttot-50, squeeze(sigfper_saveLR(1,ps,:,exrun))/Ne, 'color', sett.purple, 'linewidth', 2)
%     plot([pulse_start, pulse_start], [5,6.5], 'color', [243,161,0]/255, 'linewidth', 2)
%     ylabel('Spike rate (Hz)')
%     ylim([5,6.5])
%     xlim(xl)
%
%     subplot(413); hold on
%     plot(dthist:dthist:sett.Ttot-0.5, squeeze(sigfreq0_saveLR(1,1,:,exrun)), 'color', [0.6, 0.6,0.6], 'linewidth', 2)
%     plot(dthist:dthist:sett.Ttot-0.5, squeeze(sigfreq_saveLR(1,ps,:,exrun)), 'color', sett.blue, 'linewidth', 2)
%     plot([pulse_start, pulse_start], [70,80], 'color', [243,161,0]/255, 'linewidth', 2)
%     ylabel('Frequency (Hz)')
%     xlim(xl)
%
%     subplot(414); hold on
%     plot(dthist:dthist:sett.Ttot, squeeze(sigpower0_saveLR(1,1,:,exrun)./1000.^2), 'color', [0.6, 0.6,0.6], 'linewidth', 2)
%     plot(dthist:dthist:sett.Ttot, squeeze(sigpower_saveLR(1,ps,:,exrun)./1000.^2), 'color', sett.green, 'linewidth', 2)
%     plot([pulse_start, pulse_start], [-10e-4,10e-3], 'color', [243,161,0]/255, 'linewidth', 2)
%     ylabel('Power')
%     xlim(xl)
%     xlabel('Time (ms)')
% end

%% new measures

% figure;
% gray = [0.5, 0.5, 0.5];
%
% subplot(511); hold on
% plot(phases, h_spikes(1,Ids,1), 'Color', gray);
% title('Spikes in first period after pulse')
% ylabel('# spikes')
% ylim([0, 400])
%
%
% subplot(512); hold on
% plot(phases, h_fper(1,Ids,1), 'Color', gray);
% title('Highest concentration of new spikes');
% ylabel('# spikes')
% ylim([0, 800])
%
% subplot(513); hold on
% plot(phases, h_base(1,Ids,1), 'Color', sett.purple);
% title('Low-pass filtered STH')
% ylabel('Firing rate change (au)')
% ylim([0, 12])
%
% subplot(514); hold on
% plot(phases, h_freq(1,Ids,1), 'Color', sett.green);
% title('Frequency change')
% ylabel('Hz')
% ylim([-2, 4])
%
% subplot(515); hold on
% plot(phases, h_peak(1,Ids,1), 'Color', sett.blue);
% title('Difference in amplitude');
% ylabel('au')
% xlabel('Phase of pulse')
% ylim([25,175])
%

%% deltaHeights & deltaTimes
% avgPhase = nanmean(sort(phases,1),2);
% avgCV = mean(CVind_e, 1);
%
%
% figure('Position', [50,50, 200,600]);
% subaxis(1,1,1, 'SpacingVert',0.05,'SpacingHoriz', 0.05,'MarginTop',0.1,'MarginBottom',0.10, 'MarginRight',0.2,'MarginLeft',0.4, 'Padding',0);hold on
% scatter(ones(length(avgCV),1), avgCV, 40, colors(1:length(avgCV),:), 'filled');
% ylabel('CV_{norm}')
%
% for e = 1:esteps
%
%     % heights
%     [peaktemp_plus, timetemp_plus] = findpeaks(heights1(:,e), 'sort', 'descend', 'npeaks', 1, 'minPeakHeight', 50);
%     if ~isempty(peaktemp_plus)
%         ytemp = [heights1(timetemp_plus:end,e); heights1(1:timetemp_plus-1,e)];
%     else
%         ytemp = heights1(:,e);
%     end
%
%     [peaktemp_min, timetemp_min] = findpeaks(-heights1(:,e), 'sort', 'descend', 'npeaks', 1, 'minPeakHeight', 0);
%     if ~isempty(peaktemp_min)
%         ytemp = [heights1(timetemp_min:end,e); heights1(1:timetemp_min-1,e)];
%     else
%         ytemp = heights1(:,e);
%     end
%
%     figure(3); subplot(211); hold on
%     plot(phases(:,e), heights1(:,e), 'color', colors(e,:));
%     plot(phases(timetemp_plus,e), peaktemp_plus, 'o', 'color', colors(e,:), 'MarkerFaceColor', colors(e,:));
%     plot(phases(timetemp_min,e), -peaktemp_min, '^', 'color', colors(e,:), 'MarkerFaceColor', colors(e,:));
%
%     subplot(212); hold on;
%     plot(phases(timetemp_plus,e), peaktemp_plus, 'o', 'color', colors(e,:), 'MarkerFaceColor', colors(e,:));
%     plot(phases(timetemp_min,e), -peaktemp_min, '^', 'color', colors(e,:), 'MarkerFaceColor', colors(e,:));
%
%     % shifts
%     [peaktemp_plus, timetemp_plus] = findpeaks(shifts1(:,e), 'sort', 'descend', 'npeaks', 1, 'minPeakHeight', -2);
%     if ~isempty(peaktemp_plus)
%         ytemp = [shifts1(timetemp_plus:end,e); shifts1(1:timetemp_plus-1,e)];
%     else
%         ytemp = shifts1(:,e);
%     end
%     [peaktemp_min, timetemp_min] = findpeaks(-shifts1(:,e), 'sort', 'descend', 'npeaks', 1, 'minPeakHeight', 2);
%     if ~isempty(peaktemp_min)
%         ytemp = [shifts1(timetemp_min:end,e); shifts1(1:timetemp_min-1,e)];
%     else
%         ytemp = shifts1(:,e);
%     end
%
%     figure(4); subplot(211); hold on
%     plot(phases(:,e),shifts1(:,e), 'color', colors(e,:));
%     plot(phases(timetemp_plus,e), peaktemp_plus, 'o', 'color', colors(e,:), 'MarkerFaceColor', colors(e,:));
%     plot(phases(timetemp_min,e), -peaktemp_min, '^', 'color', colors(e,:), 'MarkerFaceColor', colors(e,:));
%
%     subplot(212); hold on;
%     plot(phases(timetemp_plus,e), peaktemp_plus, 'o', 'color', colors(e,:), 'MarkerFaceColor', colors(e,:));
%     plot(phases(timetemp_min,e), -peaktemp_min, '^', 'color', colors(e,:), 'MarkerFaceColor', colors(e,:));
%
% end
