function [psd, f] = spectraldensity(hista, Fs, reg) 

nfft = pow2(nextpow2(length(hista(:,:,1))));
psd = zeros(nfft/2+1, 1, reg);

for r = 1:reg
[psd(:,:,r), f] = pmtm(hista(:,:,r), 2, nfft, Fs);
end

end