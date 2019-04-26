% INFO and CORRELATIONS

l_amp = arrayfun(@(x) length(sett.Iistep_amp{x}), 1:length(sett.Iistep_amp));
nsteps = max(l_amp);

for ns = 1:nsteps
    
    %% Freq and Power
    [freqosc_st, ~, power_st, ~] = freq_phase(hist_e(1,(Tselection + (ns-1)*sett.Istep_dur)/dthist+1:...
        (Tselection + ns*sett.Istep_dur)/dthist,:), dthist, sett.corselection, reg,[]);
    
    sigfreq = freqosc_st(1,100:end-99,:);
    sigpower = power_st(1,100:end-99,:);
    
    %% Firing rate
    % firing rate per period
    
    sigfpertmp = nan(1,sett.Istep_dur/dthist - 199, reg);
    histtmp = hist_e(1,(Tselection+(ns-1)*sett.Istep_dur)/dthist+1:(Tselection+ns*sett.Istep_dur)/dthist,:);
    
    for r = 1:reg
        
        iw = 1;
        wnd = 0;
        avghist = NaN(1,size(histtmp,2));
        while iw+wnd <= size(histtmp,2)
            freqtmp = freqosc_st(1,iw,r);
            wnd = round((1000/freqtmp)/dthist);  
            if isnan(wnd)
                wnd = round((1000/nanmean(freqosc_st(1,:,r)))/dthist);
            end
            if iw+wnd > size(histtmp,2)
                break
            end
            avghist(1,iw) = sum(histtmp(1,iw:iw+wnd-1,r));
            iw = iw+1;
        end
        sigfpertmp(1,1:length(avghist),r) = smooth(avghist(1,1:end), round((2000/nanmean(freqosc_st(1,:,r)))/dthist));
        
%         if r == 3
%             addwin = -1;%0;
%             addsmooth = 1;%+200;
%         else
%             addwin = -1;
%             addsmooth = 1;
%         end
        
%         freqtmp = min(freqosc_st(1,:,r)); %nanmean(freqosc_st(1,:,r), 2);
%         wnd = round((1000/freqtmp)/dthist) + addwin; %-1;
%         avghist = zeros(1,size(histtmp,2)-wnd);
%         for iw = 1:size(histtmp,2)-wnd
%             avghist(1,iw) = sum(histtmp(1,iw:iw+wnd,r));
%         end
%         sigfper(1,:,r) = smooth(avghist(1,100-wnd/2:end-100+wnd/2), round((1000/freqtmp)/dthist)+addsmooth);
    end
%     
    sigfper = sigfpertmp(1,100:end-100,:);
%     figure; hold on
%     plot(freqosc_st(1,:,1), 'b')
%     plot(freqosc_st(1,:,2), 'g')
%     plot(freqosc_st(1,:,3), 'r')
%     
%     figure; hold on;
%     plot(avghist(1,:), 'r')
%     
%     figure; hold on;
%     plot(sigfper(1,:,1), 'b')
%     plot(sigfper(1,:,2), 'g')
%     plot(sigfper(1,:,3), 'r')
    
    %% normalize
    
    sigfper_norm = nan(1,sett.Istep_dur/dthist - 199, reg);
    sigfreq_norm = nan(1,sett.Istep_dur/dthist - 199, reg);
    sigpower_norm = nan(1,sett.Istep_dur/dthist - 199, reg);
    
    for r = 1:reg
        sigfreq_norm(1,:,r) = (sigfreq(1,:,r) - min(sigfreq(1,:,r))) / max(sigfreq(1,:,r)-min(sigfreq(1,:,r)));
        sigpower_norm(1,:,r) = (sigpower(1,:,r) - min(sigpower(1,:,r))) / max(sigpower(1,:,r)- min(sigpower(1,:,r)));
        sigfper_norm(1,:,r) = (sigfper(1,:,r) - min(sigfper(1,:,r))) / max(sigfper(1,:,r)- min(sigfper(1,:,r)));
    end
    
    %% interreg info transfer - corr
    for r = 1:reg-1
        reg1 = r;
        reg2 = reg;
        [rho_freq_part(ns,r), p_freq_part(ns,r)] = nancorr(squeeze(sigfreq_norm(1,:,reg1))', squeeze(sigfreq_norm(1,:,reg2))');
        [rho_power_part(ns,r), p_power_part(ns,r)] = nancorr(squeeze(sigpower_norm(1,:,reg1))', squeeze(sigpower_norm(1,:,reg2))');
        [rho_fper_part(ns,r), p_fper_part(ns,r)] = nancorr(squeeze(sigfper_norm(1,:,reg1))', squeeze(sigfper_norm(1,:,reg2))');
    end
    
end

