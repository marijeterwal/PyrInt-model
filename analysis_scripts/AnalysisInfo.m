% INFO and CORRELATIONS

if sum(sett.Iestep) + sum(sett.Iistep)
    AnalysisInfoStepWise
else
    
    % normalization of noise signal
    sig_noise = nsig(:,Tselection/dthist:end-100);
    for r = 1:reg
        sig_noise_norm(r,:) = (sig_noise(r,:) - min(sig_noise(r,:))) ./ max(sig_noise(r,:)- min(sig_noise(r,:)));
    end
    
    % preallocate
    sigfreq = freqosc(1,Tselection/dthist:end-99,:);
    sigpower = power(1,Tselection/dthist:end-99,:);
    % sigpeak = zeros(1,(Ttot-Tselection)/dthist - 99, reg);
    % sigbase = zeros(1,(Ttot-Tselection)/dthist - 99, reg);
    sigfper = zeros(1,(Ttot-Tselection)/dthist - 99, reg);
    
    Fs = 1000/dthist;
    % [b a] = butter(5, [30,90]./(Fs/2), 'stop'); % 10th order low pass filter, below 30 Hz
    % [b a] = butter(5, 200./(Fs/2), 'high'); % 10th order low pass filter, below 30 Hz
    % [b a] = butter(5, 20./(Fs/2), 'low'); % 10th order low pass filter, below 30 Hz
    % fvtool(b,a)
    
    [sigfreq_tshuf, sigpower_tshuf, sigfper_tshuf, ...
        sigfreq_pshuf, sigpower_pshuf, sigfper_pshuf] = deal(zeros(size(sigfreq)));
    
    for r = 1:reg
        
        %     % peak signal
        %     s = 0.25*(1000/avgmorfreq(l1,l2,r))/dthist;
        %     c = s*1*sqrt(2*pi);
        %     gp = (c/(s*sqrt(2*pi)))*exp(-(-50:dthist:50).^2/(2*s^2));
        %
        %     deltaT = round(100/(avgmorfreq(l1,l2,r)*dthist));
        %     [peaks, times] = findpeaks(squeeze(hist_e(1,:,r)), 'sortstr', 'descend', 'npeaks', round(Ttot/floor(1000/avgmorfreq(l1,l2,r))), 'minpeakdistance', round(500/(dthist*avgmorfreq(l1,l2,r))));
        %     [times, Id] = sort(times);
        %     peaks = peaks(Id);
        %
        %     psig = zeros(1,length(hist_e(1,:,1)));
        %     psig(times) = peaks;
        %     psig2 = conv(psig, gp, 'same');
        %     sigpeak(1,:,r) = psig2(1,Tselection/dthist:end-100);
        
        % baseline firing rate
        %     [b a] = butter(5, [avgmorfreq(l1,l2,r)-20,avgmorfreq(l1,l2,r)+50]./(Fs/2), 'stop'); % 10th order low pass filter, below 30 Hz
        %     bsig = filtfilt(b, a, hist_e(:,:,r));
        %     bsig = hist_e(:,:,r);
        %     bsig = smooth(bsig,round((1000/avgmorfreq(l1,l2,r))/dthist));
        %     sigbase(1,:,r) = bsig(Tselection/dthist:end-100);
        
        % firing rate per period
        wnd = round((1000/avgmorfreq(l1,l2,r))/dthist)-1;
        avghist = zeros(1,size(hist_e,2)-wnd);
        for iw = 1:size(hist_e,2)-wnd
            avghist(1,iw) = sum(hist_e(1,iw:iw+wnd,r));
        end
        sigfper(1,:,r) = smooth(avghist(1,Tselection/dthist-wnd/2:end-100+wnd/2), round((1000/avgmorfreq(l1,l2,r))/dthist));

        % normalize signals
        %     sigbase_norm(1,:,r) = (sigbase(1,:,r) - min(sigbase(1,:,r))) / max(sigbase(1,:,r)- min(sigbase(1,:,r)));
        sigfreq_norm(1,:,r) = (sigfreq(1,:,r) - min(sigfreq(1,:,r))) / max(sigfreq(1,:,r)- min(sigfreq(1,:,r)));
        sigpower_norm(1,:,r) = (sigpower(1,:,r) - min(sigpower(1,:,r))) / max(sigpower(1,:,r)- min(sigpower(1,:,r)));
        %     sigpeak_norm(1,:,r) = (sigpeak(1,:,r) - min(sigpeak(1,:,r))) / max(sigpeak(1,:,r)- min(sigpeak(1,:,r)));
        sigfper_norm(1,:,r) = (sigfper(1,:,r) - min(sigfper(1,:,r))) / max(sigfper(1,:,r)- min(sigfper(1,:,r)));
        
        % calculate own baseline - correlations
        [rho_freq(l1,l2,r), p_freq(l1,l2,r)] = nancorr(sig_noise(r,:)', squeeze(sigfreq(1,:,r))');
        [rho_power(l1,l2,r), p_power(l1,l2,r)] = nancorr(sig_noise(r,:)', squeeze(sigpower(1,:,r))');
        %     [rho_peak(l1,l2,r), p_peak(l1,l2,r)] = nancorr(sig_noise(r,:)', squeeze(sigpeak(1,:,r))');
        %     [rho_base(l1,l2,r), p_base(l1,l2,r)] = nancorr(sig_noise(r,:)', squeeze(sigbase(1,:,r))');
        [rho_fper(l1,l2,r), p_fper(l1,l2,r)] = nancorr(sig_noise(r,:)', squeeze(sigfper(1,:,r))');
        
        %     for l = 1:length(20:100)
        %     noisemat = sig_noise(r,:);%repmat(sig_noise(r,:), [length(powermat(:,1,1)),1]);
        %     powermat2 = squeeze(mean(powermat(5:16,Tselection/dthist:end-100,r),1));
        %         [rho_powermat(l1,l2,r), ~] = corr(noisemat(:), powermat2(:));
        %     end
        %     rho_powermat(l1,l2,r) = max(rho_temp);
        
        %     % mutual info
        windowsizeMI = 1000/avgmorfreq(l1,l2,r);
        nMI_freq(l1,l2,r) = mutualinformation(sig_noise_norm(r,:), squeeze(sigfreq_norm(1,:,r)), 0, 1, nbinsMI);
        nMI_power(l1,l2,r) = mutualinformation(sig_noise_norm(r,:), squeeze(sigpower_norm(1,:,r)), 0, 1, nbinsMI);
        %     nMI_peak(l1,l2,r) = mutualinformation(sig_noise(r,:)', squeeze(sigpeak(1,:,r))', 0, 1, nbinsMI);
        %     nMI_base(l1,l2,r) = mutualinformation(sig_noise_norm(r,:), squeeze(sigbase_norm(1,:,r)), 0, 1, nbinsMI);
        nMI_fper(l1,l2,r) = mutualinformation(sig_noise_norm(r,:), squeeze(sigfper_norm(1,:,r)), 0, 1, nbinsMI);
        
    end
    
    % noise vs signal - normalized
    for r = 1:reg-1
        reg1 = r;
        reg2 = reg;
        
        [rho_freq(l1,l2,reg+r), p_freq(l1,l2,reg+r)] = nancorr(sig_noise(reg1,:)', squeeze(sigfreq(1,:,reg2))');
        [rho_power(l1,l2,reg+r), p_power(l1,l2,reg+r)] = nancorr(sig_noise(reg1,:)', squeeze(sigpower(1,:,reg2))');
        %     [rho_peak(l1,l2,reg+r), p_peak(l1,l2,reg+r)] = nancorr(sig_noise(reg1,:)', squeeze(sigpeak(1,:,reg2))');
        %     [rho_base(l1,l2,reg+r), p_base(l1,l2,reg+r)] = nancorr(sig_noise(reg1,:)', squeeze(sigbase(1,:,reg2))');
        [rho_fper(l1,l2,reg+r), p_fper(l1,l2,reg+r)] = nancorr(sig_noise(reg1,:)', squeeze(sigfper(1,:,reg2))');
        
        %     for l = 1:length(20:100)
        %     noisemat = sig_noise(reg1,:);%repmat(sig_noise(reg1,:), [length(powermat(:,1,1)),1]);
        %     powermat2 = squeeze(mean(powermat(5:16,Tselection/dthist:end-100,reg2), 1));
        %         [rho_powermat(l1,l2,reg+r), ~] = corr(noisemat(:), powermat2(:));
        %     end
        %     rho_powermat(l1,l2,reg+r) = max(rho_temp);
        
        nMI_freq(l1,l2,reg+r) = mutualinformation(sig_noise_norm(reg1,:), squeeze(sigfreq_norm(1,:,reg2)), 0, 1, nbinsMI);
        nMI_power(l1,l2,reg+r) = mutualinformation(sig_noise_norm(reg1,:), squeeze(sigpower_norm(1,:,reg2)), 0, 1, nbinsMI);
        %     nMI_peak(l1,l2,reg+r) = mutualinformation(sig_noise_norm(reg1,:)', squeeze(sigpeak(1,:,reg2))', 0, 1, nbinsMI);
        %     nMI_base(l1,l2,reg+r) = mutualinformation(sig_noise_norm(reg1,:), squeeze(sigbase_norm(1,:,reg2)), 0, 1, nbinsMI);
        nMI_fper(l1,l2,reg+r) = mutualinformation(sig_noise_norm(reg1,:), squeeze(sigfper_norm(1,:,reg2)), 0, 1, nbinsMI);
        
    end
    
    %% between areas
    
    for nc = 1:ncomb
        reg1 = comb(nc,1);
        reg2 = comb(nc,2);
        hbins = round(1000/avgmorfreq(l1,l2,reg1)); % half a period
        
        % prep signals for MI: remove nan
        sigpn = sigpower_norm; sigpn(isnan(sigpn)) = 0;
        sigfn = sigfreq_norm; sigfn(isnan(sigfn)) = 0;
        sigfpn = sigfper_norm; sigfpn(isnan(sigfpn)) = 0;
        
        MI_power_12(l1,l2,nc) = mutualinformation(squeeze(sigpn(1,:,reg1)), squeeze(sigpn(1,:,reg2)), 1, hbins, nbinsMI);
        MI_freq_12(l1,l2,nc) = mutualinformation(squeeze(sigfn(1,:,reg1)), squeeze(sigfn(1,:,reg2)), 1, hbins, nbinsMI);
        MI_fper_12(l1,l2,nc) = mutualinformation(squeeze(sigfpn(1,:,reg1)), squeeze(sigfpn(1,:,reg2)), 1, 1, nbinsMI);
        
        
        rho_power_12(l1,l2,nc) = nancorr(squeeze(sigpower_norm(1,:,reg1))', squeeze(sigpower_norm(1,:,reg2))');
        % rho_peak_12(l1,l2,nc) = nancorr(squeeze(sigpeak_norm(1,:,reg1))', squeeze(sigpeak_norm(1,:,reg2))');
        % rho_base_12(l1,l2,nc) = nancorr(squeeze(sigbase_norm(1,:,reg1))', squeeze(sigbase_norm(1,:,reg2))');
        rho_freq_12(l1,l2,nc) = nancorr(squeeze(sigfreq_norm(1,:,reg1))', squeeze(sigfreq_norm(1,:,reg2))');
        rho_fper_12(l1,l2,nc) = nancorr(squeeze(sigfper_norm(1,:,reg1))', squeeze(sigfper_norm(1,:,reg2))'); 
        
    end
end

%% shuffled


% time-shuffled signals
sigfreq_tshuf(:,:,1) = sigfreq(:,:,1);
sigpower_tshuf(:,:,1) = sigpower(:,:,1);
sigfper_tshuf(:,:,1) = sigfper(:,:,1);

sigfreq_tshuf(:,:,2) = sigfreq(:,randperm(size(sigfreq,2)),2);
sigpower_tshuf(:,:,2) = sigpower(:,randperm(size(sigfreq,2)),2);
sigfper_tshuf(:,:,2) = sigfper(:,randperm(size(sigfreq,2)),2);

% period-shuffled signals
wnd = round((1000/avgmorfreq(l1,l2,r))/dthist)-1; % period window
nwnd = floor(size(sigfper,2)/wnd);
randwnd = randperm(nwnd);

sigfreq_pshuf = sigfreq;
sigpower_pshuf = sigpower;
sigfper_pshuf = sigfper;

for iw = 1:nwnd
    sigfreq_pshuf(1,(iw-1)*wnd+1:iw*wnd,2) = sigfreq(1, (randwnd(iw)-1)*wnd+1:randwnd(iw)*wnd,2);
    sigpower_pshuf(1,(iw-1)*wnd+1:iw*wnd,2) = sigpower(1, (randwnd(iw)-1)*wnd+1:randwnd(iw)*wnd,2);
    sigfper_pshuf(1,(iw-1)*wnd+1:iw*wnd,2) = sigfper(1, (randwnd(iw)-1)*wnd+1:randwnd(iw)*wnd,2);
end


rho_power_12t(l1,l2) = nancorr(squeeze(sigpower(1,:,1))', squeeze(sigpower_tshuf(1,:,2))');
rho_freq_12t(l1,l2) = nancorr(squeeze(sigfreq(1,:,1))', squeeze(sigfreq_tshuf(1,:,2))');
rho_fper_12t(l1,l2) = nancorr(squeeze(sigfper(1,:,1))', squeeze(sigfper_tshuf(1,:,2))'); 

rho_power_12p(l1,l2) = nancorr(squeeze(sigpower(1,:,1))', squeeze(sigpower_pshuf(1,:,2))');
rho_freq_12p(l1,l2) = nancorr(squeeze(sigfreq(1,:,1))', squeeze(sigfreq_pshuf(1,:,2))');
rho_fper_12p(l1,l2) = nancorr(squeeze(sigfper(1,:,1))', squeeze(sigfper_pshuf(1,:,2))');
        