function [spikes] = findspikes(Y, dt, threshold)
% Y is membrane voltage in format neurons x times

N = length(Y(:,1,1));
tnumb = length(Y(1,:,1));
reg = length(Y(1,1,:));

% spikeseries = zeros(1,tnumb, reg);
spikes = [0,0,0];

for r = 1:reg
for t = 2:tnumb
    for n = 1:N
        if (Y(n,t,r) >= threshold) && (Y(n,t-1,r) <= threshold)
            spikes = [spikes; n, (t-1) * dt,r]; %t-1 because at location 1, t = 0;
%             spikeseries(1,t,r) = 1;
        end
    end
end
end
spikes = spikes(2:end,:);
end