function MI = mutualinformation(sig1, sig2, removeZeros, dt, bins)
%this code takes two signals (time domain, row vectors): sig1 and sig2
%dt: the time step used for horizontal binning, in time indices
%bins: either an integer number of vertical bins or the bin edges

assert(length(sig1) == length(sig2), 'Input signals must have equal sizes');

switch nargin
    case 5
        ChangeHor = 1;
        ChangeVer = 1;
    case 4
        ChangeHor = 1;
        ChangeVer = 0;
    case 3
        ChangeHor = 0;
        ChangeVer = 0;
end

% change binsize
if ChangeHor == 1
    hbins = floor(length(sig1)/dt);
    Hsig1 = zeros(1,hbins);
    Hsig2 = zeros(1,hbins);
    for i = 1:hbins
        if i*dt > length(sig1(1,:))
            iend = length(sig1(1,:));
        else
            iend = i*dt;
        end
        Hsig1(1,i) = sum(sig1(1,(i-1)* dt +1: iend));
        Hsig2(1,i) = sum(sig2(1,(i-1)* dt +1: iend));
    end
else
    Hsig1 = Sig1;
    Hsig2 = Sig2;
end


% figure; hold on
% % plot(sig1, 'b');
% plot(Hsig1,'r')
% plot(Hsig2,'b')

% determine bin edges
if ChangeVer == 1
    if length(bins) == 1
        binmax1 = max(Hsig1);
        binmin1 = min(Hsig1);
        binmax2 = max(Hsig2);
        binmin2 = min(Hsig2);
        
        binEdges1 = linspace(binmin1, binmax1, bins+1);
        binEdges2 = linspace(binmin2, binmax2, bins+1);
        binSeries = 1:bins;
        nbins = bins;
    else
        binEdges = bins;
        binEdges2 = bins;
        nbins = length(bins)-1;
        binSeries = 1:bins;
    end
    
    Dsig1 = zeros(1,hbins);
    Dsig2 = zeros(1,hbins);
    for d = 1:nbins
        Dsig1(1,Hsig1 >= binEdges1(d) & Hsig1 <= binEdges1(d+1)) = d;
        Dsig2(1,Hsig2 >= binEdges2(d) & Hsig2 <= binEdges2(d+1)) = d;
    end
else
    Dsig1 = Hsig1;
    Dsig2 = Hsig2;
    binmax = max([Hsig1,Hsig2]);
    binmin = min([Hsig1,Hsig2]);
    binSeries = binmin:binmax;
end

if removeZeros
    Id = ~ismember(1:numel(Dsig1), find(Dsig1 == 0 & Dsig2 == 0));
    Dsig1d = Dsig1(Id);
    Dsig2d = Dsig2(Id);
else
    Dsig1d = Dsig1;
    Dsig2d = Dsig2;
end

% figure; hold on
% plot(Dsig1d, 'b');
% plot(Dsig2d, 'r');

Px = hist1D(Dsig1d, binSeries);
Py = hist1D(Dsig2d, binSeries);
Pxy = hist2D(Dsig1d, Dsig2d, binSeries);

% disp(Px)
% disp(Py)
% disp(Pxy)
% figure; imagesc(Pxy); set(gca, 'ydir', 'normal');

MI = mutinfo(Px, Py, Pxy);
end

function N = hist1D(sig, bins)
%hist for digitized data

nbins = length(bins);

N = zeros(nbins, 1);

for i = 1:nbins
    N(i) = sum(sig==bins(i));
end

N = N / length(sig);
% figure;
% bar(bins, N);
% set(gca, 'YDir', 'normal');

end

function N = hist2D(sig1, sig2, bins)
%hist for digitized data

nbins = length(bins);

N = zeros(nbins, nbins);

for i = 1:nbins
    tempsig = sig2(sig1 == bins(i));
    for j = 1:nbins
        N(j,i) = sum(tempsig==bins(j));
    end
end

N = N /length(sig1);

% figure;
% imagesc(bins, bins, N);
% set(gca, 'YDir', 'normal');

end

function MI = mutinfo(Px, Py, Pxy)

assert(length(Px) == length(Py) && length(Py) == length(Pxy(1,:)) && length(Py) == length(Pxy(:,1)), 'Input must have equal sizes');

nbins = length(Px);
MI = 0;

for n = 1:nbins %y
    for k = 1:nbins %x
        if Pxy(n,k) ~= 0
            MI = MI + Pxy(n,k)*log2(Pxy(n,k)/(Px(k) *Py(n)));
        end
    end
end
end


