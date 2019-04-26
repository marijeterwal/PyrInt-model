function [delay] = drivedelay(corei, dt)

[heights, ~] = findpeaks(corei);%,
if isempty(heights)
    delay = NaN;
else
    [~, times] = findpeaks(corei, 'minpeakheight', 0.1);% minpeakheigt is required to prevent from finding small bumps in the correlation below 0
    half = ceil(length(corei)/2);
    times = times - half;
    if isempty(times); delay = NaN;
    else [~, I] = sort(abs(times));
        times = times(I);
        delay = times(1)*dt;
    end
end
end