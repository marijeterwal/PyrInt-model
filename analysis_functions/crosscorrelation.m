function [cor, cov] = crosscorrelation(hist1, hist2, max_delay)

% assert(numel(hist1) == numel(hist2), 'Inputs for xcor must be of equal size.');

tnumb = length(hist1(1,:));
if isempty(max_delay); max_delay = tnumb;
elseif max_delay + 1 > tnumb; max_delay = tnumb;
else max_delay = max_delay +1; 
end

covplus = zeros(1,tnumb-1); % TODO should be max_delay
corplus = zeros(1,tnumb-1);
covmin = zeros(1,tnumb-1);
cormin = zeros(1,tnumb-1);

avg1 = mean(hist1(1,:));
avg2 = mean(hist2(1,:));
std1 = std(hist1(1,:));
std2 = std(hist2(1,:));

for i = 1 : max_delay
    covmin(i) = sum((hist2(1, 1:(end-(i-1)))- avg2) .* (hist1(1,i:end)- avg1));
    cormin(i) = sum((hist2(1, 1:(end-(i-1))) - avg2) .* (hist1(1,i:end) - avg1)) / (std1 *std2);
    
    covplus(i) = sum((hist1(1, 1:(end-(i-1)))- avg1) .* (hist2(1,i:end)- avg2));
    corplus(i) = sum((hist1(1, 1:(end-(i-1))) - avg1) .* (hist2(1,i:end) - avg2)) / (std1 *std2);
end

cov = [fliplr(covmin(1:end)), covplus(2:end)];
cor = [fliplr(cormin(1:end)), corplus(2:end)];

cov = cov/length(cov);
cor = cor/length(cor);

end