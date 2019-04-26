
% noise traces

tsett = sett;

if exist('csett', 'var') == 1
    names = fieldnames(csett);
    for n = 1:length(names)
        tsett.(names{n}) = csett.(names{n});
    end
end

s = RandStream('mcg16807','Seed',tsett.noiseseed);
RandStream.setGlobalStream(s);

% noiseind = (sum(tsett.Iinoise) + sum(tsett.Ienoise));
%Brown
tau = tsett.tau; %in ms
alpha = 1-dt/tau;
afil = [1 -alpha];
bfil = 1;

nsigb = zeros(reg, tnumb);
for r = 1:reg
    xtb = tsett.noiseamp * 0.5*randn(1, tnumb);
    xfb = filter(bfil,afil,xtb)*sqrt(1-alpha);
    nsigb(r,:) = xfb;
end
nsigb(reg,:) = tsett.noisecor*nsigb(1,:) + (1-tsett.noisecor)*nsigb(reg,:);

step = dthist/dt;
nsig = zeros(reg, round(length(nsigb(1,:))/step));
for t = 1:round(length(nsigb(1,:))/step)
    nsig(:,t) = (0.5+mean(nsigb(:,(t-1)*step+1:t*step),2)); % 300*
end