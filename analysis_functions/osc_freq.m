function [freq] = osc_freq(cor, dt, reg)
freq = zeros(reg,1);

for r = 1:reg
[heights, ~] = findpeaks(cor(:,:,r), 'sortstr', 'descend');%,
if isempty(heights)
    freq(r) = NaN;
else
%     threshold = heights(1)/2;
threshold = 0;
    [~, times] = findpeaks(cor(:,:,r), 'minpeakheight', threshold); %findpeaks does not consider 0 , 'sortstr', 'descend'
    if isempty(times); freq = NaN;
    else freq(r) = 1000 / (times(1)*dt);
    end
end
end
end