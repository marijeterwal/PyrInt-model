function [histt, histtime, smoothhist] = histograms(Y, N, Ttot, dt, dthist, bandpass_low, bandpass_high, reg)
% Y is spikesmatrix

steps = round((Ttot/dt +1) / (dthist / dt));
histt = zeros(1, steps, reg);
histtime = 0.5*dthist:dthist:(Ttot-0.5*dthist);
smoothhist = zeros(1, steps, reg);

Fs = 1000/(dthist);

%lowpass
fcutsl = [bandpass_high (bandpass_high + 10)];
magsl = [1 0];
devsl = [0.01 0.05];
[nl,Wnl,betal,ftypel] = kaiserord(fcutsl,magsl,devsl,Fs);
hhl = fir1(nl,Wnl,ftypel,kaiser(nl+1,betal),'noscale');

%highpass
lowband_low = bandpass_low-10;
if lowband_low < 0; lowband_low = 0; end
fcutsh = [lowband_low bandpass_low];
magsh = [0 1];
devsh = [0.05 0.01];
[nh,Wnh,betah,ftypeh] = kaiserord(fcutsh,magsh,devsh,Fs);
hhh = fir1(nh,Wnh,ftypeh,kaiser(nh+1,betah),'noscale');

hhtot = conv(hhl, hhh, 'same');

for r = 1:reg
    data = Y(Y(:,3) == r,2);
    histt(1,:,r) = (1000/(N*dthist)) * hist(data, histtime);
    smoothhist(:,:,r) = filtfilt(hhtot,1,histt(1,:,r)); %zero-phase filtering
end

smoothhist(smoothhist < 0) = 0;

end