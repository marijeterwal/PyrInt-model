function [peakheight, peaktime] = determinecorpeak(cor, dt)
% 1 = autocor; 2 = crosscor
% slope is expressed in s^-1

nr = 2;
time = zeros(nr,1);
height = zeros(nr,1);

cordata = zeros(1,floor(length(cor(1,:,1))/2)+1, nr);
cordata(:,:,1) = fliplr(cor(1: floor(length(cor(1,:,1))/2)+1)); %delay < 0
cordata(:,:,2) = cor(floor(length(cor(1,:,1))/2)+1 : end); % delay > 0

for r = 1:nr
    
    peaksnumb = 5;
%     threshold = 0;
    [threshold, ttime] = findpeaks(cor, 'sortstr', 'descend', 'npeaks', 1);
    if ttime == floor(length(cor(1,:,1))/2)+1; 
        heights = threshold; 
        times = 0;
    else
        [heights, times] = findpeaks(cordata(:,:,r), 'minpeakheight', threshold/2, 'npeaks', peaksnumb); %findpeaks does not consider 0 , 'sortstr', 'descend'
    end
    
    if isempty(times); times(r) = NaN; heights(r) = NaN;
    else
        if r == 1; times = - times * dt; end
        if r == 2; times = times *dt; end
    end
        time(r) = times(1);
        height(r) = heights(1);
end

if isnan(time(1)) && isnan(time(2))
    peaktime = NaN;
    peakheight = NaN;
else
    if abs(time(1)) == abs(time(2))
        if time(1) == 0; 
            id = 1;
        else
        id = height == max(height(1), height(2));
        end
    else
        id = abs(time) == min(abs(time(1)), abs(time(2)));
    end
    peaktime = time(id)*dt;
    peakheight = height(id);
end
end
