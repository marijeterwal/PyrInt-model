function rate = firingrate(spikes,N,reg)
% spikes is a matrix of number of spikes x 2. Column 1 registers the neuron
% number, column 2 the spiketime
rate = zeros(reg,1);


if N == 1
    for r = 1:reg
        data_ind = spikes(spikes(:,3) == r, :);
        rate(r) = 1000*length(data_ind)/(max(data_ind(:,2))-min(data_ind(:,2)));
    end
else
for r = 1:reg
    data_ind = spikes(spikes(:,3) == r, :);
    
    data = zeros(N,1);
    rtemp = zeros(N,1);
    for m = 1:N
        spikessorted = data_ind(data_ind(:,1) == m,2); %spiketimes of neuron m
        if length(spikessorted) < 2 ; continue; end
        data(m) = mean(diff(spikessorted)); %mean interspike interval of neuron m
        rtemp(m) = 1/data(m);
    end
    rtemptot = sum(rtemp);
    rate(r) = 1/N * rtemptot *1000; %from ms^-1 to s^-1, data ~= 0 is used to exclude neurons with less than 2 spikes
end
end
end