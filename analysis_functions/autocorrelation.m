function [cor, cov] = autocorrelation(hista, max_delay, reg)

tnumb = length(hista(1,:));
if isempty(max_delay); max_delay = tnumb; 
elseif max_delay + 1 > tnumb; max_delay = tnumb;
else max_delay = max_delay +1; 
end

cov = zeros(1, max_delay, reg);
cor = zeros(1, max_delay, reg);
histb = hista;

for r = 1:reg
avg = mean(histb(1,:,r));
stdev = std(histb(1,:,r));
covt = zeros(1, max_delay);
cort = zeros(1, max_delay);
for i = 1 : max_delay;
    covt(i) = sum((histb(1, 1:(end-(i-1)), r)- avg) .* (histb(1, i:end, r)- avg));
    cort(i) = sum((histb(1, 1:(end-(i-1)), r) - avg) .* (histb(1, i:end, r) - avg)) / (stdev^2);
end
cov(:,:,r) = covt/length(histb);
cor(:,:,r) = cort/length(histb);
end

end