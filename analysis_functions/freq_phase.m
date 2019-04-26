function [freq, phase, power, powermat] = freq_phase(sig, dthist, corselection, reg, reffreq)

if isempty(reffreq)
%     [coree, ~] = autocorrelation(sig, corselection/dthist, reg);
    freqtemp = avg_freq(sig, dthist, reg);
else
    freqtemp = reffreq;
end

phase = zeros(1,length(sig(:,:,1)), reg);
freq = zeros(1,length(sig(:,:,1))-1, reg);
power = zeros(1,length(sig(:,:,1))-1, reg);

for r = 1:reg
    
%     freqtemp(r)
    
    %% determine phase at reffreq
    scale = 1000/(freqtemp(r)*dthist); %scale is based on oscillation frequency and sampling frequency: Fs/Fosc
    coef = cwt(sig(:,:,r), scale-1:scale+1, 'cmor1-1');
    phase(1,:,r) = - (angle(mean(coef,1)) / (2*pi)); %reversed and shifted definition of phase...
    phase(phase < 0) = 1 + phase(phase < 0);
    
    %% phase --> frequency
    freqFromPhase = diff(unwrap(phase(:,:,r)*2*pi)) /(2*pi*dthist/1000);
    
%         figure;
%         plot(freqFromPhase)
    
    %% exclude abnormal behaviours
    excl = find(abs(diff(diff(unwrap(phase(1,:,r)*2*pi)))) >= 5*median(abs(diff(diff(unwrap(phase(1,:,r)*2*pi))))));
    exclall = excl;
    for t = 1:ceil(500/freqtemp(r)/dthist) % exlude period around anomaly
        exclall = [exclall, excl-t, excl+t;];
    end

    freq(:,:,r) = smooth(freqFromPhase, 500/freqtemp(r)/dthist); %freqFromPhase;
    
    if ~isempty(exclall)
        exclall = sort(exclall);
        exclall = exclall(exclall>0 & exclall<=length(sig(:,:,1))-1);
        freq(:,exclall,r) = NaN;
    end

    freq(1,freq(1,:,r)<=0,r) = NaN;
    
    %% power at frequency
    
    % scales of interest
    scales = unique(1000./(freq(1,:,r)*dthist));
    scales(end+1) = 0;
    tmp = 1000./(freq(1,:,r)*dthist); % scale at each time
    tmp(isnan(tmp)) = 0;
    
    % scales --> scaleid
    scaletrace = arrayfun(@(x) find(tmp(x) == scales), 1:size(freq,2));
    scaletrace(scaletrace == length(scales)) = NaN;
    
    % wavelets
    coefmat = abs(cwt(sig(:,:,r), scales(1:end-1), 'cmor1-1'));
    powermat(:,:,r) = coefmat;
    freqscales = scal2frq(scales(1:end-1), 'cmor1-1', dthist/1000); % (1000/dthist)./ scale;
    
    % scaleid --> power
    idx = sub2ind(size(coefmat),scaletrace,1:length(tmp));
    power(1,~isnan(idx),r) = coefmat(idx(~isnan(idx)));
    power(1,isnan(idx),r) = NaN;
    
    %% plot
%     figure;
%     subplot(4,1,1:3);
%     hold on;
%     pcolor((1:size(coefmat,2))*dthist ,freqscales, abs(coefmat));
%     plot((1:size(coefmat,2)-1)*dthist,freq(:,:,r), 'r')
%     shading flat
%     
%     subplot(4,1,4)
%     plot((1:size(coefmat,2)-1)*dthist, power(1,:,r))
    
end
end


function [freq] = avg_freq(sig, dt, reg)
freq = zeros(reg,1);

for r = 1:reg
    cortmp = xcorr(sig(1,:,r));
    cor = smooth(cortmp(1,round(length(cortmp)/2):end),5);
    
    [heights, ~] = findpeaks(cor, 'sortstr', 'descend');%,
    if isempty(heights)
        freq(r) = NaN;
    else
        threshold = 0;
        [~, times] = findpeaks(cor, 'minpeakheight', threshold, 'minpeakdistance', 10/dt, 'npeaks', 10); %findpeaks does not consider 0 , 'sortstr', 'descend'
        if isempty(times); freq = NaN;
        else
            difftimes = diff(times);
            freq(r) = 1000 / (mean(difftimes)*dt);
        end
    end
end
end



% compute freq from phase and smooth
%     s = 2; %ms
%     gp = (dthist/s)*(1/(sqrt(2*pi)))*exp(-(-5*s:dthist:5*s).^2/(2*s^2));
%     freq(:,:,r) = conv(freqFromPhase, gp, 'same'); % the old way
%     power(1,:,r) = conv(abs(coef), gp, 'same');

% z = hilbert(hist_e(1,:,2));
% dS = 5;
% Fs = (1000/0.5)/dS
% az = angle(z);
% figure; plot(az)
%
% instfreq = Fs/(2*pi)*diff(unwrap(2*az(1:dS:end)));%Fs/(2*pi)*diff(unwrap(angle(z)));
% S = 10;
% figure; plot(smooth(instfreq,S));
%
% figure
% spectrogram(hist_e(1,:,1),200,180,256,2000,'yaxis')
% shading interp
%
% figure
% spectrogram(hist_e(1,:,2),200,180,256,2000,'yaxis')
% shading interp

%     scale = 1000/(freqtemp(1)*dthist);
%     if 10 < scale && scale < 200 % ??
%     scalemat = scale-30:scale+300;
%     freqscales = scal2frq(scalemat, 'cmor1-1', dthist/1000); % (1000/dthist)./ scale;
%     coefmat = cwt(sig(:,:,r), scalemat, 'cmor1-1');
%     powermat(:,:,r) = abs(coefmat);




% function freq  = freqinoscillation(sig, dthist, deltafreq, reg)
% %sig is a histogram with dimensions 1 x Ttot/dthist x reg
% % flipid = 1 if phase needs to be assessed from data before timepoint
%
% freq = zeros(1,length(sig(:,:,1)), reg);
%
% for r = 1:reg
%     scale = 10:200; %10-200 Hz;
%     freqscales = scal2frq(scale, 'cmor1-1', dthist/1000); % (1000/dthist)./ scale;
%
%     coef = cwt(sig(:,:,r), scale, 'cmor1-1');
%
%     [peak,time] = findpeaks(mean(abs(coef(freqscales>15 & freqscales<150,:)),2), 'sortstr', 'descend', 'npeaks', 1, 'minpeakdistance',5);
%     if numel(peak) == 0
%         freq(1,:,r) = NaN(1,length(sig(1,:,r)),1);
%     else
%         freqest = freqscales(time + find(freqscales<150,1) - 1);
%
%     minID = find(freqscales<freqest+deltafreq,1);
%     for l = 1:length(sig(1,:,r))
%         [heights, peaks] = findpeaks(abs(coef(freqscales>freqest-deltafreq & freqscales<freqest+deltafreq,l)), 'sort', 'descend');
%         if numel(heights) == 0
% %             peakcoef(l,r) = NaN;
% %             peakscale(l) = NaN;
%             freq(1,l,r) = NaN;
%         else
%             freq(1,l,r) = freqscales(minID + peaks(1) - 1);
% %         else
% %             for i = 1:10;
% %                 if numel(heights)<i; break; end
% %                 if freqscales(peaks(i)) < freqest-30 || freqscales(peaks(i)) > freqest+30
% %                     peakcoef(l,r) = NaN;
% %                     peakscale(l) = NaN;
% %                 else
% %                     peakcoef(l,r) = heights(i);
% %                     peakscale(l) = scale(peaks(i));
% %                     break
% %                 end
% %             end
%         end
%     end
%     end
% %     freqcoefs(:,:,r) = abs(coef);
% %
% %     freq(1,:,r) = scal2frq(peakscale, 'cmor1-1', dthist/1000); %(1000/dthist)./ peakscale(:);
%     peakcoef = 0;
%     freqcoefs = 0;
%
%     %        figure;
%     %         subplot(411); plot(sig(:,:,r));
%     %         subplot(412); plot(abs(coef(1,:)));
%     %         subplot(413); imagesc(1:length(sig(1,:,r)), scale, abs(coef));
%     %         subplot(414); plot(freq(1,:,r)); ylim([0 100]);
%
% end
% end

