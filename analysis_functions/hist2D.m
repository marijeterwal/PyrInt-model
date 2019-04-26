function N = hist2D(sig1, sig2, dt, bins)

assert(length(sig1) == length(sig2), 'Input signals must have equal sizes');

% determine bin edges
if length(bins) == 1
    binmax1 = max(sig1);
    binmin1 = min(sig1);
    binmax2 = max(sig2);
    binmin2 = min(sig2);

    binEdges1 = linspace(binmin1, binmax1, bins+1);
    binEdges2 = linspace(binmin2, binmax2, bins+1);
    nbins = bins;
else
    binEdges1 = bins;
    binEdges2 = bins;
    nbins = length(bins)-1;
end

%digitize signals
hbins = round(length(sig1)/dt);
Hsig1 = zeros(1,hbins);
Hsig2 = zeros(1,hbins);
for i = 1:hbins
    if i*dt > length(sig1(1,:)) 
        iend = length(sig1(1,:));
    else 
        iend = i*dt;
    end
    Hsig1(1,i) = mean(sig1(1,(i-1)* dt +1: iend));
    Hsig2(1,i) = mean(sig2(1,(i-1)* dt +1: iend));
end


Dsig1 = zeros(1,hbins);
Dsig2 = zeros(1,hbins);
for d = 1:nbins
    Dsig1(1,Hsig1 >= binEdges1(d) & Hsig1 <= binEdges1(d+1)) = d;
    Dsig2(1,Hsig2 >= binEdges2(d) & Hsig2 <= binEdges2(d+1)) = d;
end

binsSeries = 1:nbins;

N = zeros(nbins, nbins);

for i = 1:nbins
    tempsig = Dsig2(Dsig1 == binsSeries(i));
    for j = 1:nbins
            N(j,i) = sum(tempsig==binsSeries(j));
    end
end

% figure;
% imagesc(bins, bins, N);
% set(gca, 'YDir', 'normal');

end