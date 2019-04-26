

%% Coherence in bins

if reg>1
    
    nthres = 10;
    
    cohbins = linspace(0,1,21);
    rho_freq_col = nan(length(cohbins)-1, ncomb);
    rho_peak_col = nan(length(cohbins)-1, ncomb);
    rho_base_col = nan(length(cohbins)-1, ncomb);
    
    diff_freq = nan(length(cohbins)-1, ncomb);
    diff_peak = nan(length(cohbins)-1, ncomb);
    diff_base = nan(length(cohbins)-1, ncomb);
    
    diff_freq_std = nan(length(cohbins)-1, ncomb);
    diff_peak_std = nan(length(cohbins)-1, ncomb);
    diff_base_std = nan(length(cohbins)-1, ncomb);
    
    for d = 1:length(cohbins)-1;
        for nc = 1:ncomb
            sig1 = sigfreqtot(1,cohtot(nc,:)>cohbins(d) & cohtot(nc,:)<=cohbins(d+1), comb(nc,1));
            sig2 = sigfreqtot(1,cohtot(nc,:)>cohbins(d) & cohtot(nc,:)<=cohbins(d+1), comb(nc,2));
            if isempty(sig1) || length(sig1) < nthres
                continue
            end
            [rho_freq_col(d,nc), p_freq_col(d,nc)] = corr(sig1', sig2');
            diff_freq(d,nc) = mean((sig1 - sig2)./sig2);
            diff_freq_std(d,nc) = std((sig1 - sig2)./sig2);
            
            sig1 = sigpeaktot(1,cohtot(nc,:)>cohbins(d) & cohtot(nc,:)<=cohbins(d+1), comb(nc,1));
            sig2 = sigpeaktot(1,cohtot(nc,:)>cohbins(d) & cohtot(nc,:)<=cohbins(d+1), comb(nc,2));
            [rho_peak_col(d,nc), p_peak_col(d,nc)] = corr(sig1', sig2');
            diff_peak(d,nc) = mean((sig1 - sig2)./sig2);
            diff_peak_std(d,nc) = std((sig1 - sig2)./sig2);
            
            sig1 = sigbasetot(1,cohtot(nc,:)>cohbins(d) & cohtot(nc,:)<=cohbins(d+1), comb(nc,1));
            sig2 = sigbasetot(1,cohtot(nc,:)>cohbins(d) & cohtot(nc,:)<=cohbins(d+1), comb(nc,2));
            [rho_base_col(d,nc), p_base_col(d,nc)] = corr(sig1', sig2');
            diff_base(d,nc) = mean((sig1 - sig2)./sig2);
            diff_base_std(d,nc) = std((sig1 - sig2)./sig2);
        end
    end
    
    xax = 0+diff(cohbins(1:2))/2:diff(cohbins(1:2)):1-diff(cohbins(1:2))/2;
    
    rho_freq_sig = rho_freq_col; rho_freq_sig(p_freq_col > 0.05) = NaN;
    rho_peak_sig = rho_peak_col; rho_peak_sig(p_freq_col > 0.05) = NaN;
    rho_base_sig = rho_base_col; rho_base_sig(p_freq_col > 0.05) = NaN;
    
    figure;
    subplot(131); hold on
    plot(xax, rho_freq_sig(:,2), 'color', sett.purple);%, 30, sett.purple, 'filled')
    plot(xax, rho_freq_sig(:,1), 'color', sett.green);%, 30, sett.green, 'filled')
    
    subplot(132); hold on
    plot(xax, rho_peak_sig(:,2), 'color', sett.purple);%, 30, sett.purple, 'filled')
    plot(xax, rho_peak_sig(:,1), 'color', sett.green);%, 30, sett.green, 'filled')
    
    subplot(133); hold on
    plot(xax, rho_base_sig(:,2), 'color', sett.purple);%, 30, sett.purple, 'filled')
    plot(xax, rho_base_sig(:,1), 'color', sett.green);%, 30, sett.green, 'filled')
    
    
    
    figure;
    subplot(131); hold on
    plot(xax, mean(rho_freq_col(:,1:2), 2), 'k');%, 30, sett.purple, 'filled')
    xlim([0,1]); ylim([-1, 1]); title('LFP frequency')
    ylabel('Correlation coefficient');
    
    subplot(132); hold on
    plot(xax, mean(rho_peak_col(:,1:2), 2), 'k');%, 30, sett.purple, 'filled')
    xlim([0,1]); ylim([-1, 1]); title('LFP amplitude')
    xlabel('Coherence');
    
    subplot(133); hold on
    plot(xax, mean(rho_base_col(:,1:2), 2), 'k');%, 30, sett.purple, 'filled')
    xlim([0,1]); ylim([-1, 1]); title('Firing rate')
    
    
    figure; hold on
    subplot(131); hold on
    errorbar(xax, diff_freq(:,2), diff_freq_std(:,2), 'color', sett.purple)
    errorbar(xax, diff_freq(:,1), diff_freq_std(:,1), 'color', sett.green)
    xlim([0,1]);
    
    subplot(132); hold on
    errorbar(xax, diff_peak(:,2), diff_peak_std(:,2), 'color', sett.purple)
    errorbar(xax, diff_peak(:,1), diff_peak_std(:,1), 'color', sett.green)
    xlim([0,1]);
    
    subplot(133); hold on
    errorbar(xax, diff_base(:,2), diff_base_std(:,2), 'color', sett.purple)
    errorbar(xax, diff_base(:,1), diff_base_std(:,1), 'color', sett.green)
    xlim([0,1]);
    
    %% Phase in bins
    
    phasedata = phasetot;
    phasedata(1,cohtot < 0.8) = NaN;
    
    phasebins = linspace(0,1,41);
    delays = 0;%(0:windowstepsize:60)/windowstepsize;
    
    rho_freq_col = nan(length(phasebins)-1, length(delays), ncomb);
    rho_peak_col = nan(length(phasebins)-1, length(delays), ncomb);
    rho_base_col = nan(length(phasebins)-1, length(delays), ncomb);
    
    p_freq_col = nan(length(phasebins)-1, length(delays), ncomb);
    p_peak_col = nan(length(phasebins)-1, length(delays), ncomb);
    p_base_col = nan(length(phasebins)-1, length(delays), ncomb);
    
    diff_freq = nan(length(phasebins)-1, length(delays), ncomb);
    diff_peak = nan(length(phasebins)-1, length(delays), ncomb);
    diff_base = nan(length(phasebins)-1, length(delays), ncomb);
    
    diff_freq_std = nan(length(phasebins)-1, length(delays), ncomb);
    diff_peak_std = nan(length(phasebins)-1, length(delays), ncomb);
    diff_base_std = nan(length(phasebins)-1, length(delays), ncomb);
    
    % for d = delays;
    d = delays(1);
    sigfreq_shortd = sigfreqtot;%sigfreq_short(1,d+1:end,:);
    sigpeak_shortd = sigpeaktot;%sigpeak_short(1,d+1:end,:);
    sigbase_shortd = sigbasetot;%sigbase_short(1,d+1:end,:);
    phasedatad = phasedata;%phasedata(1,1:end-d,:);
    
    for p = 1:length(phasebins)-1;
        for nc = 1:ncomb
            sig1 = sigfreq_shortd(1,phasedatad(1,:,nc)>phasebins(p) & phasedatad(1,:,nc)<=phasebins(p+1), comb(nc,1));
            sig2 = sigfreq_shortd(1,phasedatad(1,:,nc)>phasebins(p) & phasedatad(1,:,nc)<=phasebins(p+1), comb(nc,2));
            if isempty(sig1) || length(sig1) < 5
                continue
            end
            diff_freq(p,d+1,nc) = mean((sig1 - sig2)./sig2);
            diff_freq_std(p,d+1,nc) = std((sig1 - sig2)./sig2);
            [rho_freq_col(p,d+1,nc),p_freq_col(p,d+1,nc)] = corr(sig1', sig2');
            
            sig1 = sigpeak_shortd(1,phasedatad(1,:,nc)>phasebins(p) & phasedatad(1,:,nc)<=phasebins(p+1), comb(nc,1));
            sig2 = sigpeak_shortd(1,phasedatad(1,:,nc)>phasebins(p) & phasedatad(1,:,nc)<=phasebins(p+1), comb(nc,2));
            diff_peak(p,d+1,nc) = mean((sig1 - sig2)./sig2);
            diff_peak_std(p,d+1,nc) = std((sig1 - sig2)./sig2);
            [rho_peak_col(p,d+1,nc),p_peak_col(p,d+1,nc)] = corr(sig1', sig2');
            
            sig1 = sigbase_shortd(1,phasedatad(1,:,nc)>phasebins(p) & phasedatad(1,:,nc)<=phasebins(p+1), comb(nc,1));
            sig2 = sigbase_shortd(1,phasedatad(1,:,nc)>phasebins(p) & phasedatad(1,:,nc)<=phasebins(p+1), comb(nc,2));
            diff_base(p,d+1,nc) = mean((sig1 - sig2)./sig2);
            diff_base_std(p,d+1,nc) = std((sig1 - sig2)./sig2);
            [rho_base_col(p,d+1,nc),p_base_col(p,d+1,nc)] = corr(sig1', sig2');
        end
    end
    % end
    
    xax = 0+diff(phasebins(1:2))/2:diff(phasebins(1:2)):1-diff(phasebins(1:2))/2;
    d = 1;
    
    rho_freq_sig = rho_freq_col; rho_freq_sig(p_freq_col > 0.05) = NaN;
    rho_peak_sig = rho_peak_col; rho_peak_sig(p_freq_col > 0.05) = NaN;
    rho_base_sig = rho_base_col; rho_base_sig(p_freq_col > 0.05) = NaN;
    
    figure;
    subplot(131); hold on
    plot(xax, rho_freq_col(:,d,2), 'color', sett.purple);%, 30, sett.purple, 'filled')
    plot(xax, rho_freq_col(:,d,1), 'color', sett.green);%, 30, sett.green, 'filled')
    
    subplot(132); hold on
    plot(xax, rho_peak_col(:,d,2), 'color', sett.purple);%, 30, sett.purple, 'filled')
    plot(xax, rho_peak_col(:,d,1), 'color', sett.green);%, 30, sett.green, 'filled')
    
    subplot(133); hold on
    plot(xax, rho_base_col(:,d,2), 'color', sett.purple);%, 30, sett.purple, 'filled')
    plot(xax, rho_base_col(:,d,1), 'color', sett.green);%, 30, sett.green, 'filled')
    
    figure;
    subplot(131); hold on
    scatter(xax, rho_freq_col(:,d,2),30, sett.purple, 'filled')
    scatter(xax, rho_freq_col(:,d,1), 30, sett.green, 'filled')
    
    subplot(132); hold on
    scatter(xax, rho_peak_col(:,d,2), 30, sett.purple, 'filled')
    scatter(xax, rho_peak_col(:,d,1), 30, sett.green, 'filled')
    
    subplot(133); hold on
    scatter(xax, rho_base_col(:,d,2), 30, sett.purple, 'filled')
    scatter(xax, rho_base_col(:,d,1), 30, sett.green, 'filled')
    
    
    figure;
    subplot(131); hold on
    plot(xax, mean(rho_freq_col(:,1:2), 2), 'k');%, 30, sett.purple, 'filled')
    xlim([0,1]); ylim([-0.8, 0.8]); title('LFP frequency')
    ylabel('Correlation coefficient');
    
    subplot(132); hold on
    plot(xax, mean(rho_peak_col(:,1:2), 2), 'k');%, 30, sett.purple, 'filled')
    xlim([0,1]); ylim([-0.8, 0.8]); title('LFP amplitude')
    xlabel('Phase difference');
    
    subplot(133); hold on
    plot(xax, mean(rho_base_col(:,1:2), 2), 'k');%, 30, sett.purple, 'filled')
    xlim([0,1]); ylim([-0.8, 0.8]); title('Firing rate')
    
    
    figure;
    subplot(131); scatter(xax, mean(diff_freq(:,d,1:2),3), 30, 'k', 'filled'); xlim([0,1])
    subplot(132); scatter(xax, mean(diff_peak(:,d,1:2),3), 30, 'k', 'filled'); xlim([0,1])
    subplot(133); scatter(xax, mean(diff_base(:,d,1:2),3), 30, 'k', 'filled'); xlim([0,1])
    
    figure;
    subplot(131); plot(xax, mean(diff_freq(:,d,1:2),3),'k')
    subplot(132); plot(xax, mean(diff_peak(:,d,1:2),3),'k')
    subplot(133); plot(xax, mean(diff_base(:,d,1:2),3),'k')
    
    
    figure; hold on
    subplot(231); hold on; xlim([0,1]);
    errorbar(xax, diff_freq(:,d,2), diff_freq_std(:,d,2), 'color', sett.purple)
    subplot(234); hold on; xlim([0,1]);
    errorbar(xax, diff_freq(:,d,1), diff_freq_std(:,d,1), 'color', sett.green)
    
    subplot(232); hold on; xlim([0,1]);
    errorbar(xax, diff_peak(:,d,2), diff_peak_std(:,2), 'color', sett.purple)
    subplot(235); hold on; xlim([0,1]);
    errorbar(xax, diff_peak(:,d,1), diff_peak_std(:,1), 'color', sett.green)
    
    subplot(233); hold on; xlim([0,1]);
    errorbar(xax, diff_base(:,d,2), diff_base_std(:,d,2), 'color', sett.purple)
    subplot(236); hold on; xlim([0,1]);
    errorbar(xax, diff_base(:,d,1), diff_base_std(:,d,1), 'color', sett.green)
    
    figure; hold on
    subplot(131); hold on; xlim([0,1]);
    scatter(xax, diff_freq(:,d,2), 30, sett.purple, 'filled')
    scatter(xax, diff_freq(:,d,1), 30, sett.green, 'filled')
    
    subplot(132); hold on; xlim([0,1]);
    scatter(xax, diff_peak(:,d,2), 30, sett.purple, 'filled')
    scatter(xax, diff_peak(:,d,1), 30, sett.green, 'filled')
    
    subplot(133); hold on; xlim([0,1]);
    scatter(xax, diff_base(:,d,2), 30, sett.purple, 'filled')
    scatter(xax, diff_base(:,d,1), 30, sett.green, 'filled')
    
end