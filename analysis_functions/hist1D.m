function N = hist1D(sig, dt, bins)
%hist for digitized data

% determine bin edges
if length(bins) == 1
    binmax = max(sig);
    binmin = min(sig);

    binEdges = linspace(binmin, binmax, bins+1);
    nbins = bins;
else
    binEdges = bins;
    nbins = length(bins)-1;
end

%digitize signals
hbins = round(length(sig)/dt);
Hsig = zeros(1,hbins);
for i = 1:hbins
    if i*dt > length(sig(1,:)) 
        iend = length(sig(1,:));
    else 
        iend = i*dt;
    end
    Hsig(1,i) = mean(sig(1,(i-1)* dt +1: iend));
end

Dsig = zeros(1,hbins);
for d = 1:nbins
    Dsig(1, Hsig > binEdges(d) & Hsig < binEdges(d+1)) = d;
end

N = zeros(nbins, 1);

for i = 1:nbins
    N(i) = sum(Dsig==i);
end

% figure;
% bar(bins, N);
% set(gca, 'YDir', 'normal');

end