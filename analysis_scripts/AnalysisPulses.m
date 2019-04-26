% analyze effects of pulses

pulse_start = csett.Ipulse_start;

phaseosc1 = phaseosc(1,Tselection/dthist:end-100,:);

phaseOfPulse(l1,l2) = phaseosc0_t(pulse_start/dthist);

% phase shifting
shift = zeros(size(phaseosc1));

t1 = (pulse_start - Tselection)/dthist+1;
t2 = (pulse_start - Tselection + 100)/dthist +1;

for r = 1:reg
    
    phaseshiftseries = (unwrap(phaseosc1(:,:,r)*2*pi) - unwrap(phaseosc0(:,:,r)*2*pi))/ (2*pi);
    local_phaseshift(l1,l2,r) = wrapTo2Pi(mean(phaseshiftseries(1,t1:t2))*2*pi)/(2*pi);
    local_phaseshift_full(l1,l2,r) = wrapTo2Pi(mean(phaseshiftseries(1,t1:end-100))*2*pi)/(2*pi);
    local_phaseshift_std(l1,l2,r) = wrapTo2Pi(std(phaseshiftseries(1,t1:end-100))*2*pi)/(2*pi);
end

if reg > 1
    
    numbercomb = nchoosek(reg,2);
    comb = flipud(combnk(1:reg,2));
    
    for nc = 1:numbercomb
        reg1 = comb(nc,1);
        reg2 = comb(nc,2);
        
        phaseshiftseries = (unwrap(phaseosc1(:,:,reg1)*2*pi) - unwrap(phaseosc1(:,:,reg2)*2*pi))/ (2*pi);
        reg_phaseshift(l1,l2,nc) = wrapTo2Pi(mean(phaseshiftseries(1,t1:t2))*2*pi)/(2*pi);
        reg_phaseshift_full(l1,l2,nc) = wrapTo2Pi(mean(phaseshiftseries(1,t1:end-100))*2*pi)/(2*pi);
        reg_phaseshift_std(l1,l2,nc) = wrapTo2Pi(std(phaseshiftseries(1,t1:end-100)) *2*pi)/(2*pi);
    end
    
end

for r = 1:reg
    
    startPlot = pulse_start - 50;
    endPlot = startPlot + 200;
    
    %% peaks
    
    % Determine smoothhiste peak heights and times: npeaks +1 peaks
    pulseTime = pulse_start; % ms
    % smoothhist or hist?
    data = smoothhist_e(1,pulse_start/dthist-1:round((pulse_start + (npeaks+1.1)*1000/avgmorfreq(l1,l2,1))/dthist),r);
    [heights_t, times_t] = findpeaks(data, 'npeaks', npeaks+1, 'minPeakDistance', round((0.6*1000/avgmorfreq(l1,l2,1))/dthist)); % 'sortstr', 'descend',
    times_t = (times_t-1)*dthist + pulse_start -dthist; % ms
    [temptimes, Id] = sort(times_t);
    times = [temptimes, NaN(1,npeaks+1 - length(times_t))];
    heights = [heights_t, NaN(1,npeaks+1 - length(times_t))];
    
    % Compare peaks with no-pulse trace; reduce to npeaks peaks
    if pulseTime >= times0(1,1)
        times0_temp = times0(r,2:end);
        heights0_temp = heights0(r,2:end);
        times_temp = times(2:end);
        heights_temp = heights(2:end);
    else
        times0_temp = times0(r,1:npeaks);
        heights0_temp = heights0(r,1:npeaks);
        times_temp = times(1:npeaks);
        heights_temp = heights(1:npeaks);
    end
    deltaTimes(l1,l2,r,:) = times_temp - times0_temp; % positive is later
    deltaHeights(l1,l2,r,:) = heights_temp - heights0_temp;
    
    % more analysis
    %         tempsig = sigpeak(1,(startPlot)/dthist:(endPlot)/dthist,r)- sigpeak0(1,(startPlot)/dthist:(endPlot)/dthist,r);
    %         [~,time] = findpeaks(abs(tempsig), 'sortstr', 'descend', 'npeaks', 1); %
    %         if isempty(time);
    %             h_peak(l1,l2,r) = NaN;
    %         else
    %             h_peak(l1,l2,r) = tempsig(time);
    %         end
    if r == 2 && times0_temp(1) <= pulse_start + sett.delay_reg_e + 2.5
        h_peakheight0(l1,l2,r) = heights0_temp(2);
        h_peakheight1(l1,l2,r) = heights_temp(2);
        
        h_peakheight(l1,l2,r) = deltaHeights(l1,l2,r,2);
    else
        h_peakheight0(l1,l2,r) = heights0_temp(1);
        h_peakheight1(l1,l2,r) = heights_temp(1);
        
        h_peakheight(l1,l2,r) = deltaHeights(l1,l2,r,1);
    end
    
    if r == 2 && times0_temp(1) <= pulse_start + sett.delay_reg_e + 2.5
        h_peaktime0(l1,l2,r) = times0_temp(2);
        h_peaktime1(l1,l2,r) = times_temp(2);
        
        h_peaktime(l1,l2,r) = deltaTimes(l1,l2,r,2);
    else
        h_peaktime0(l1,l2,r) = times0_temp(1);
        h_peaktime1(l1,l2,r) = times_temp(1);
        
        h_peaktime(l1,l2,r) = deltaTimes(l1,l2,r,1);
    end
    
    if h_peaktime(l1,l2,r) > (1000/avgmorfreq(l1,l2,r))/2
        h_peaktime(l1,l2,r) = h_peaktime(l1,l2,r) - (1000/avgmorfreq(l1,l2,r));
    end
    
    
    %         s = 1*(1000/avgmorfreq(l1,l2,r))/dthist;
    %         c = s*1*sqrt(2*pi);
    %         gp = (c/(s*sqrt(2*pi)))*exp(-(-100:dthist:100).^2/(2*s^2));
    %
    %         deltaT = round(100/(avgmorfreq(l1,l2,r)*dthist));
    %         [peaks, times] = findpeaks(squeeze(hist_e(1,:,r)), 'sortstr', 'descend', 'npeaks', round(Ttot/floor(1000/avgmorfreq(l1,l2,r))), 'minpeakdistance', round(500/(dthist*avgmorfreq(l1,l2,r))));
    %         [times, Id] = sort(times);
    %         peaks = peaks(Id);
    %
    %         psig = zeros(1,length(hist_e(1,:,1)));
    %         psig(times) = peaks;
    %         psig2 = conv(psig, gp, 'same');
    %         sigpeak(1,:,r) = psig2;
    
    [peaks, times] = findpeaks(squeeze(smoothhist_e(1,(pulse_start)/dthist:...
        (pulse_start)/dthist+round((1000/avgmorfreq(l1,l2,r))/dthist),r)), 'npeaks', 1, 'sortstr', 'descend');
    sigpeak(r) = peaks;
    
    h_peak(l1,l2,r) = sigpeak(r) - sigpeak0(r);
    
    %% firing rates
    
    % sliding window firing rate
    wnd = round((1000/avgmorfreq(l1,l2,r))/dthist)-1;
    avghist = zeros(1,size(hist_e,2)-wnd);
    for iw = 1:size(hist_e,2)-wnd
        avghist(1,iw) = sum(hist_e(1,iw:iw+wnd,r));
    end
    sigfper(1,:,r) = smooth(avghist(1,Tselection/dthist-wnd/2:end-100+wnd/2), round((1000/avgmorfreq(l1,l2,r))/dthist));
    
    [height, ~] = findpeaks(sigfper(1,round((startPlot-Tselection)/dthist):round((endPlot-Tselection)/dthist),r) - sigfper0(1,round((startPlot-Tselection)/dthist):round((endPlot-Tselection)/dthist),r), 'sortstr', 'descend', 'npeaks', 1);
    if isempty(height);
        h_fper(l1,l2,r) = NaN;
    else
        h_fper(l1,l2,r) = height;
    end
    
    % difference in spiking in one period after pulse
    h_spikes0(l1,l2,r) = sum(hist_e0(1,round(pulse_start/dthist):round(pulse_start/dthist+1*(1000/avgmorfreq0(1))/dthist),r), 2);
    h_spikes1(l1,l2,r) = sum(hist_e(1,round(pulse_start/dthist):round(pulse_start/dthist+1*(1000/avgmorfreq0(1))/dthist),r), 2);
    
    h_spikes(l1,l2,r) = sum(hist_e(1,round(pulse_start/dthist):round(pulse_start/dthist+1*(1000/avgmorfreq0(1))/dthist),r)-hist_e0(1,round(pulse_start/dthist):round(pulse_start/dthist+1*(1000/avgmorfreq0(1))/dthist),r), 2);
    
    
    Fs = 1000/dthist;
    [z p k] = butter(5, 20./(Fs/2), 'low'); % 5th order low pass filter, below 20 Hz
    [sos,g] = zp2sos(z,p,k); % Convert to 2nd order sections form
    
    bsig = filtfilt(sos, g, hist_e(1,:,r));
    sigbase(1,:,r) = bsig;
    [height,~] = findpeaks(sigbase(1,(startPlot)/dthist:(endPlot)/dthist,r)- sigbase0(1,(startPlot)/dthist:(endPlot)/dthist,r), 'sortstr', 'descend', 'npeaks', 1); %
    if isempty(height);
        h_base(l1,l2,r) = NaN;
    else
        h_base(l1,l2,r) = height;
    end
    
    %% frequency
    Tper = 1000./avgmorfreq(l1,l2,r);
    [~,times] = findpeaks(power(1,round((pulse_start)/dthist):round((pulse_start+2*Tper)/dthist),r)- power0(1,round((pulse_start)/dthist):round((pulse_start+2*Tper)/dthist),r), 'sortstr', 'descend', 'npeaks', 1);
    if isempty(times)
        [~,times] = findpeaks(abs(power(1,round((pulse_start)/dthist):round((pulse_start+2*Tper)/dthist),r)- power0(1,round((pulse_start)/dthist):round((pulse_start+2*Tper)/dthist),r)), 'sortstr', 'descend', 'npeaks', 1);
    end
    if isempty(times)
        h_power(l1,l2,r) = 0;%NaN;
    else 
        h_power(l1,l2,r) = power(1, round(pulse_start/dthist) + times -1,r) - power0(1, round(pulse_start/dthist) + times -1,r);
    end
    
    %         spikesselection_e2 = spikes_e((spikes_e(:,2,:) >= startPlot) & (spikes_e(:,2,:) < endPlot), :);
    %         [CVshort, ~] = coefficientofvariation_regions(spikesselection_e2, avgmorfreq(l1,l2,:), endPlot-startPlot, reg);
    %
    %         if CVshort(r) <= 0.056
    %             h_freq(l1,l2,r) = NaN;
    %         else
    [height, ~] = findpeaks(freqosc(1,round((startPlot-Tper)/dthist):round((endPlot-Tper)/dthist),r) - freqosc0(1,round((startPlot-Tper)/dthist):round((endPlot-Tper)/dthist),r), 'sortstr', 'descend', 'npeaks', 1); %
    if isempty(height);
        h_freq(l1,l2,r) = NaN;
    else
        h_freq(l1,l2,r) = height;
%         h_freq(l1,l2,r) = height;
    end
    %         end
    
    %% synchrony
    
    spikesselection_pulse = spikesselection_e(spikesselection_e(:,2) > pulse_start &...
        spikesselection_e(:,2) < pulse_start + (1000/avgmorfreq0(r)) &...
        spikesselection_e(:,3) == r,:);
    %     phaseosc_pulse = phaseosc(1,round(pulse_start/dthist):round(pulse_start/dthist+1*(1000/avgmorfreq0(1))/dthist),r);
    sigsync = phaselocking(sett, spikesselection_pulse, phaseosc(1,:,r), [r,1]);
    
    h_sync(l1,l2,r) = sigsync - sigsync0(r);
    
    if r == 1;
       sigfreq_save(l1,l2,:) = freqosc(1,:,1);
       sigpower_save(l1,l2,:) = power(1,:,1);
       sigfper_save(l1,l2,:) = sigfper(1,:,1);
       histe_save(l1,l2,:) = hist_e(1,:,1);
    end
    
end

deltaTimes_First_Input(l1,l2,:) = deltaTimes(l1,l2,:,1);
deltaHeights_First_Input(l1,l2,:) = deltaHeights(l1,l2,:,1);
