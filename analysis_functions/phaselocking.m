
function phaselock = phaselocking(sett, spikes, reference, nc)
% Pairwise Phase Consistency
% see Vinck et al., Journal of Computational Neuroscience, 2012
% nc is a n by 2 matrix:
% - the first column indicates the regions the spikes are taken from
% - the second column indicates the regions the LFPs are taken from
% this facilitates comparing spike from different areas to the same LFP, or
% vice versa.
% if nc is empty, reg is used and area_spikes = area_lfp;

% assumed that there are no trials, or every neuron is a trial

if isempty(nc)
    nc = [(1:sett.Nregions)', (1:sett.Nregions)'];
end

phaselock = zeros(size(nc,1),1);

if size(spikes,1) > 150000
    warning('Spikes is large, PPC will take a lot of time')
end

for r = 1:size(nc,1) % comparisons

        spikeTimes = round(spikes(spikes(:,3) == nc(r,1),2)/sett.dthist); % time IDs of spikes
        if isempty(spikeTimes)
            warning('No spikes found - no PPC analysis')
        end
        phaseNr = reference(1,spikeTimes,nc(r,2))*2*pi; % phases  
        Ns = length(spikeTimes); % number of spikes

        U1 = cos(phaseNr);
        U2 = sin(phaseNr);
        U = [U1',U2'];

        Pl_n = zeros(Ns,1);
        for ns = 1:Ns % spikes
            Pl_n = Pl_n + U * [U1(ns),U2(ns)]';
        end
    phaselock(r) = (sum(Pl_n(:))-Ns)/(Ns*(Ns-1));
end

% for r = 1:nc
%     Pl_n = zeros(N,1);
%     for nn = 1:N
%         spikesNr = round(spikes(spikes(:,1) == nn & spikes(:,3) == r,2)/sett.dthist);
%         phaseNr = reference(1,spikesNr,r);
%         Pl_n(nn) = abs(sum(exp(1i*phaseNr*2*pi))/N);
%     end
%     phaselock(r) = mean(Pl_n);
% end

end