% reference run for pulses

% no pulse
fprintf('Calculating/loading data: l1 = %i of %i; l2 = %i of %i;\n', l1, l1steps, 0,l2steps);

if loaddata == 1;
    if runnumber <= 652
        if runnumber == 651 && l1 > 1
            data = open([sett.savelocdata, 'e', num2str(1),'_i', num2str(0),'.mat']);
        else
            data = open([sett.savelocdata, 'e', num2str(l1),'_i', num2str(0),'.mat']);
        end
    else
        data = open([sett.savelocdata, 'l1', num2str(l1),'_l2', num2str(0),'.mat']); % run 660
    end
    spikes_i0 = data.spikes_i0; spikes_e0 = data.spikes_e0;
else
    csett.Iipulse = zeros(1,reg);
    csett.Iepulse = zeros(1,reg);
    [spikes_i0, spikes_e0] = calcEInetworkRK4(sett, para, csett);
    %     [spikes_i0, spikes_e0] = calcEInetworkRK4(sett, para, Iinew, Ienew, zeros(Ni*reg,length(ipulsevct(1,:))), zeros(Ne*reg,length(epulsevct(1,:))), 1);
    if savedata == true; save([sett.savelocdata, 'l1', num2str(l1),'_l2', num2str(0), '.mat'], 'spikes_i0', 'spikes_e0'); end
end

[hist_i0, histtime_i0, smoothhist_i0] = histograms(spikes_i0, Ni, Ttot, dt, dthist, sett.cutoff_freq_low, sett.cutoff_freq_high, reg);
[hist_e0, histtime_e0, smoothhist_e0] = histograms(spikes_e0, Ne, Ttot, dt, dthist, sett.cutoff_freq_low, sett.cutoff_freq_high, reg);

[coree0, ~] = autocorrelation(smoothhist_e0(1,Tselection/dthist:end,:), corselection/dthist, reg);
freq_e0 = osc_freq(coree0, dthist, reg);

[freqosc0, phaseosc0_t, power0] = freq_phase(hist_e0, dthist, sett.corselection, reg, []);
avgmorfreq0 = nanmean(freqosc0(1,Tselection/dthist:end-100,:), 2);
phaseosc0 = phaseosc0_t(1,Tselection/dthist:end-100,:);

if length(sett.Ipulse_start) == 1
    
    pulse_start = sett.Ipulse_start;
    
    times0 = zeros(reg, npeaks+1);
    heights0 = zeros(reg, npeaks+1);
    for r = 1:reg
        data = smoothhist_e0(1,sett.Ipulse_start/dthist-1:round((sett.Ipulse_start + (npeaks+1.1)*1000/avgmorfreq0(1))/dthist),r);
        [heights0_t, times0_t] = findpeaks(data, 'npeaks', npeaks+1, 'minPeakDistance', round((0.6*1000/avgmorfreq0(1))/dthist));
        times0_t = (times0_t-1)*dthist + sett.Ipulse_start -dthist; % ms
        [temptimes, Id] = sort(times0_t);
        times0(r,:) = [temptimes, NaN(1,npeaks+1 - length(times0_t))];
        heights0(r,:) = [heights0_t(Id), NaN(1,npeaks+1 - length(times0_t))];
        %
%                 Tper0 = 1000./ avgmorfreq0;
        %         selvec = [Tselection/dthist:round((Ttot- Tper0(1))/dthist); round((Tselection+Tper0(1))/dthist):round((Ttot)/dthist)]';
        %         for a = 1:length(selvec(:,1))
        %             fper0(1,a,r) = sum(hist_e0(1,selvec(a,1):selvec(a,2),r), 2);
        %         end
        
        wnd = round((1000/avgmorfreq0(r))/dthist)-1;
        avghist = zeros(1,size(hist_e0,2)-wnd);
        for iw = 1:size(hist_e0,2)-wnd
            avghist(1,iw) = sum(hist_e0(1,iw:iw+wnd,r));
        end
        sigfper0(1,:,r) = smooth(avghist(1,Tselection/dthist-wnd/2:end-100+wnd/2), round((1000/avgmorfreq0(r))/dthist));
        
        Fs = 1000/dthist;
        [z p k] = butter(5, 20./(Fs/2), 'low'); % 10th order low pass filter, below 30 Hz
        [sos,g] = zp2sos(z,p,k); % Convert to 2nd order sections form
        
        %         s = 1*(1000/avgmorfreq0(r))/dthist;
        %         c = s*1*sqrt(2*pi);
        %         gp = (c/(s*sqrt(2*pi)))*exp(-(-100:dthist:100).^2/(2*s^2));
        %
        %         deltaT = round(100/(avgmorfreq0(r)*dthist));
        %         [peaks, times] = findpeaks(squeeze(hist_e0(1,:,r)), 'sortstr', 'descend', 'npeaks', round(Ttot/floor(1000/avgmorfreq0(r))), 'minpeakdistance', round(500/(dthist*avgmorfreq0(r))));
        %         [times, Id] = sort(times);
        %         peaks = peaks(Id);
        %
        %         psig = zeros(1,length(hist_e0(1,:,1)));
        %         psig(times) = peaks;
        %         psig2 = conv(psig, gp, 'same');
        %         sigpeak0(1,:,r) = psig2;
        
        [peaks, times] = findpeaks(squeeze(smoothhist_e0(1,...
            (pulse_start)/dthist:(pulse_start)/dthist+round((1000/avgmorfreq0(r))/dthist),r)),...
            'npeaks', 1, 'sortstr', 'descend');
        sigpeak0(r) = peaks;
        
        bsig = filtfilt(sos, g, hist_e0(:,:,r));
        sigbase0(1,:,r) = bsig;
        
        spikesselection_pulse = spikes_e0(spikes_e0(:,2) > pulse_start &...
            spikes_e0(:,2) < pulse_start + (1000/avgmorfreq0(r)) &...
            spikes_e0(:,3) == r,:);
        %     phaseosc_pulse = phaseosc0_t(1,round((pulse_start)/dthist):round(pulse_start/dthist+1*(1000/avgmorfreq0(1))/dthist),r);
        sigsync0(r) = phaselocking(sett, spikesselection_pulse, phaseosc0_t(1,:,r), [r,1]);
        
        if r == 1;
            sigfreq0_save(l1,l2,:) = freqosc0(1,:,1);
            sigpower0_save(l1,l2,:) = power0(1,:,1);
            sigfper0_save(l1,l2,:) = sigfper0(1,:,1);
            histe0_save(l1,l2,:) = hist_e0(1,:,1);
        end
        
    end
    sigfreq0 = freqosc0(1,Tselection/dthist:end-99,:);
else
    pulse_start = csett.Ipulse_start;
end