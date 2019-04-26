
% plot freq, sync and phase for two areas

if variableFreqAxis
    
    % account for differences in x-axis between runs
    df = avgmordata(:,1,1,:) - avgmordata(:,1,2,:);
    [sortdf,Id] = sort(df(:));
    cId = Id(~isnan(sortdf));
    [cI,cJ,cK,cL] = ind2sub(size(df),Id(~isnan(sortdf)));
    
    nf = round((length(cI)/length(runnumbers)));
    
    dfavg = zeros(nf, ncomb);
    cohavg = zeros(nf, l2steps, ncomb);
    dfreqavg = zeros(nf, l2steps, ncomb);
    phasediffavg = zeros(nf, l2steps, ncomb);
    meanPle_12avg = zeros(nf, l2steps);
    Phaseleavg = zeros(nf, l2steps);
    
    for x = 1:nf
        tmpf = zeros(1,ncomb);
        tmpcoh = zeros(1,l2steps,ncomb);
        tmpdf = zeros(1,l2steps);
        tmpp = zeros(1,l2steps,ncomb);
        tmpPPC12 = zeros(1,l2steps);
        tmpPPCi12 = zeros(1,l2steps);
        tmpPPC2 = zeros(1,l2steps);
        tmpPPCi2 = zeros(1,l2steps);
        count = 0;
        for y = 1:length(runnumbers)
            if (x-1)*length(runnumbers)+y > length(cI)
                break
            end
            tmpf = tmpf + df(cI((x-1)*length(runnumbers)+y),...
                :,:, cL((x-1)*length(runnumbers)+y));
            tmpcoh = tmpcoh + cohLR(cI((x-1)*length(runnumbers)+y),...
                :,:, cL((x-1)*length(runnumbers)+y));
            tmpdf = tmpdf + (avgmordata(cI((x-1)*length(runnumbers)+y),...
                :,2, cL((x-1)*length(runnumbers)+y)) - ...
                avgmordata(cI((x-1)*length(runnumbers)+y),...
                :,1, cL((x-1)*length(runnumbers)+y))) ;
            tmpp = tmpp + phasediffLR(cI((x-1)*length(runnumbers)+y),...
                :,:, cL((x-1)*length(runnumbers)+y));
            tmpPPC12 = tmpPPC12 + Phaselocke12LR(cI((x-1)*length(runnumbers)+y),...
                :,:, cL((x-1)*length(runnumbers)+y));
            tmpPPCi12 = tmpPPCi12 + Phaselocki12LR(cI((x-1)*length(runnumbers)+y),...
                :,:, cL((x-1)*length(runnumbers)+y));
            tmpPPC2 = tmpPPC2 + PhaselockeLR(cI((x-1)*length(runnumbers)+y),...
                :,2, cL((x-1)*length(runnumbers)+y));
            tmpPPCi2 = tmpPPCi2 + PhaselockiLR(cI((x-1)*length(runnumbers)+y),...
                :,2, cL((x-1)*length(runnumbers)+y));
            count = count + 1;
        end
        dfavg(x,:) = tmpf/count;
        cohavg(x,:) = tmpcoh/count;
        dfreqavg(x,:) = tmpdf/count;
        phasediffavg(x,:) = tmpp/count;
        meanPle_12avg(x,:) = tmpPPC12/count;
        meanPli_12avg(x,:) = tmpPPCi12/count;
        Phaseleavg(x,:) = tmpPPC2/count;
        Phaseliavg(x,:) = tmpPPCi2/count;
    end
    
    xaxis = dfavg; %avgmorfreq(:,1,1) - avgmorfreq(:,1,2);
    yaxis = sett.g_syn_ee_r(1,1):sett.g_syn_r_step:sett.g_syn_r_max;
    
    
    %% figures corrected for frequency axis
    
    figure; pcolor(xaxis, yaxis, cohavg'); axis xy;
    shading flat
    colorbar; caxis([0, 1]);
    colormap(parula)
    xlabel('\Delta frequency')
    ylabel('Total synaptic conductance (\muS/cm^{2})')
    title('Coherence')
    
    figure; pcolor(xaxis, yaxis, dfreqavg'); axis xy;
    shading flat
    colorbar; caxis([-30, 30]);
    colormap(parula)
    xlabel('\Delta frequency')
    ylabel('Total synaptic conductance (\muS/cm^{2})')
    title('\Delta freq')
    
    phasedata = phasediffavg; phasedata(cohavg<0.9) = NaN; %phasedata(phasediff_std>0.1) = NaN;%
    
    figure; pcolor(xaxis, yaxis, (phasedata*2)'); axis xy;
    shading flat
    xlabel('\Delta frequency')
    ylabel('Total synaptic conductance (\muS/cm^{2})')
    colorbar; caxis([0,2]);
    colormap(parula)
    title('Phase difference')
    
    figure;
    colormap(parula)
    subplot(121); hold on
    pcolor(xaxis, yaxis, meanPle_12avg'); axis xy; axis tight
    contour(xaxis, yaxis, round(meanPle_12avg'*100)/100, [0.25,0.50], 'w')
    % contour(xaxis, yaxis, meanPle_12avg', [0.25,0.5], 'w')
    shading flat
    title('Phase locking to area 1')
    ylabel('Total synaptic conductance (\muS/cm^{2})')
    xlabel('\Delta frequency')
    colorbar; caxis([0,1]);
    subplot(122); hold on
    pcolor(xaxis, yaxis, Phaseleavg'); axis xy; axis tight
    contour(xaxis, yaxis, Phaseleavg', [0.25,0.5], 'w')
    shading flat
    title('Phase locking to area 2')
    xlabel('\Delta frequency')
    colorbar; caxis([0,1]);
    
end

%% figure for one run

if ~variableFreqAxis
    
    xaxis = avgmorfreq(:,1,1);% - avgmorfreq(:,1,2);
    xlab = 'Frequency circuit 1';% '\Delta frequency'
    yaxis = sett.g_syn_ee_r(2,1):sett.g_syn_r_step:sett.g_syn_r_max;
    
    % figure; imagesc(xaxis, yaxis, peakcor'); axis xy; title('Xcorr peak')
    figure; pcolor(xaxis, yaxis, coh'); axis xy;
    shading flat
    colorbar; caxis([0, 1]);
    colormap(parula)
    xlabel(xlab)
    ylabel('Total synaptic conductance (\muS/cm^{2})')
    title('Coherence')
    
    figure; pcolor(xaxis, yaxis, cohi'); axis xy;
    shading flat
    colorbar; caxis([0, 1]);
    colormap(parula)
    xlabel(xlab)
    ylabel('Total synaptic conductance (\muS/cm^{2})')
    title('Coherence - interneurons')
    
    figure; pcolor(xaxis, yaxis, (avgmorfreq(:,:,2) - avgmorfreq(:,:,1))'); axis xy;
    shading flat
    colorbar; caxis([-30, 30]);
    colormap(parula)
    xlabel(xlab)
    ylabel('Total synaptic conductance (\muS/cm^{2})')
    title('\Delta freq')
    
    figure; pcolor(xaxis, yaxis, (avgmorfreqi(:,:,2) - avgmorfreqi(:,:,1))'); axis xy;
    shading flat
    colorbar; caxis([-30, 30]);
    colormap(parula)
    xlabel(xlab)
    ylabel('Total synaptic conductance (\muS/cm^{2})')
    title('\Delta freq - interneurons')
    
    % figure; pcolor(xaxis, yaxis, CVe(:,:,2)'); axis xy;
    % shading flat
    % colorbar; caxis([0.05;0.20])
    % colormap(parula)
    % xlabel('\Delta frequency')
    % ylabel('Total synaptic conductance (\muS/cm^{2})')
    % title('Resonance')
    
    % % delays & phases
    % figure; imagesc(xaxis, yaxis, delay'); axis xy;
    % title('Delay between areas')
    
    phasedata = phasediff; %phasedata(coh<0.9) = NaN; %phasedata(phasediff_std>0.1) = NaN;%
    
    figure; pcolor(xaxis, yaxis, (phasedata*2)'); axis xy;
    shading flat
    xlabel(xlab)
    ylabel('Total synaptic conductance (\muS/cm^{2})')
    colorbar; caxis([0,2]);
    colormap(parula)
    title('Phase difference')
    
    phasedata = phasediffi; phasedata(cohi<0.9) = NaN; %phasedata(phasediff_std>0.1) = NaN;%
    
    figure; pcolor(xaxis, yaxis, (phasedata*2)'); axis xy;
    shading flat
    xlabel(xlab)
    ylabel('Total synaptic conductance (\muS/cm^{2})')
    colorbar; caxis([0,2]);
    colormap(parula)
    title('Phase difference - interneurons')
    
    figure; pcolor(xaxis, yaxis, (1000./avgmorfreq(:,:,1) - delayei(:,:,1))'); axis xy;
    shading flat
    xlabel(xlab)
    ylabel('Total synaptic conductance (\muS/cm^{2})')
    colorbar; caxis([0,20]);
    colormap(parula)
    title('Delay')
    
    delaydat = delayei;
    
    figure; pcolor(xaxis, yaxis, delaydat(:,:,1)'); axis xy;
    shading flat
    xlabel(xlab)
    ylabel('Total synaptic conductance (\muS/cm^{2})')
    colorbar; caxis([-3,5]);
    colormap(parula)
    title('Delay')
    
    figure; pcolor(xaxis, yaxis, delaydat(:,:,2)'); axis xy;
    shading flat
    xlabel(xlab)
    ylabel('Total synaptic conductance (\muS/cm^{2})')
    colorbar; caxis([-3,5]);
    colormap(parula)
    title('Delay')
    
    figure; pcolor(xaxis, yaxis, (1000./avgmorfreq(:,:,2) - delayei(:,:,2))'); axis xy;
    shading flat
    xlabel(xlab)
    ylabel('Total synaptic conductance (\muS/cm^{2})')
    colorbar; caxis([0,20]);
    colormap(parula)
    title('Delay')
    
    figure;
    colormap(parula)
    subplot(121); hold on;
    pcolor(xaxis, yaxis, Phaseli(:,:,1)'); axis xy; axis tight
    shading flat
    title('Phase locking I within area 1')
    ylabel('Total synaptic conductance (\muS/cm^{2})')
    xlabel(xlab)
    colorbar; caxis([0,1]);
    subplot(122);
    pcolor(xaxis, yaxis, Phasele(:,:,1)'); axis xy;
    shading flat
    title('Phase locking E within area 1')
    xlabel('\Delta frequency')
    colorbar; caxis([0,1]);
    
    
    figure;
    colormap(parula)
    subplot(121); hold on;
    pcolor(xaxis, yaxis, meanPli_12'); axis xy; axis tight
    shading flat
    title('Phase locking to area 1')
    ylabel('Total synaptic conductance (\muS/cm^{2})')
    xlabel(xlab)
    colorbar; caxis([0,1]);
    subplot(122);
    pcolor(xaxis, yaxis, Phaseli(:,:,2)'); axis xy;
    shading flat
    title('Phase locking to area 2')
    xlabel('\Delta frequency')
    colorbar; caxis([0,1]);
    
    
    figure;
    colormap(parula)
    subplot(121); hold on
    pcolor(xaxis, yaxis, meanPle_12'); axis xy; axis tight
    contour(xaxis, yaxis, (round(meanPle_12*100)/100)', [0.25,0.5], 'w')
    shading flat
    title('Phase locking to area 1')
    ylabel('Total synaptic conductance (\muS/cm^{2})')
    xlabel(xlab)
    colorbar; caxis([0,1]);
    subplot(122); hold on
    pcolor(xaxis, yaxis, Phasele(:,:,2)'); axis xy; axis tight
    contour(xaxis, yaxis, Phasele(:,:,2)', [0.25,0.5], 'w')
    shading flat
    title('Phase locking to area 2')
    xlabel('\Delta frequency')
    colorbar; caxis([0,1]);
    
    % figure;
    % colormap(parula)
    % subplot(121); hold on
    % contour(xaxis, yaxis, meanPle_12', [0.25,0.5], 'b')
    % title('Phase locking to area 1')
    % ylabel('Total synaptic conductance (\muS/cm^{2})')
    % xlabel('\Delta frequency')
    % colorbar; caxis([0,1]);
    % subplot(122);
    % contour(xaxis, yaxis, Phasele(:,:,2)', [0.25,0.5], 'b')
    % title('Phase locking to area 2')
    % xlabel('\Delta frequency')
    % colorbar; caxis([0,1]);
    
    
    %% standard deviations
    
    if length(runnumbers) > 1
        
        figure; pcolor(xaxis, yaxis, phasediff_std'); axis xy;
        shading flat;
        xlabel(xlab)
        ylabel('Total synaptic conductance (\muS/cm^{2})')
        colorbar; caxis([0,0.1]);
        colormap(parula)
        title('Phase difference - standard deviation')
        
        figure; pcolor(xaxis, yaxis, coh_std'); axis xy;
        shading flat;
        xlabel(xlab)
        ylabel('Total synaptic conductance (\muS/cm^{2})')
        colorbar; caxis([0,0.2]);
        colormap(parula)
        title('Coherence - standard deviation')
        
        figure;
        colormap(parula)
        subplot(121); hold on
        pcolor(xaxis, yaxis, stdPle_12'); axis xy; axis tight
        shading flat
        title('Phase locking to area 1 - standard deviation')
        ylabel('Total synaptic conductance (\muS/cm^{2})')
        xlabel(xlab)
        colorbar; caxis([0,0.2]);
        subplot(122); hold on
        pcolor(xaxis, yaxis, Ple_std(:,:,2)'); axis xy; axis tight
        shading flat
        title('Phase locking to area 2 - standard deviation')
        xlabel('\Delta frequency')
        colorbar; caxis([0,0.2]);
        
    else
        
        cdata = ones(size(coh)); cdata(coh<0.9) = 0;
        
        figure;
        contour(xaxis, yaxis, cdata', [0.5]); axis xy;
        xlabel('\Delta frequency')
        ylabel('Total synaptic conductance (\muS/cm^{2})')
        colorbar; caxis([0,2]);
        colormap(parula)
        title('Coherence outline')
    end
    
end