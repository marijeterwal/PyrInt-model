function [CV_ind, CV_all] = coefficientofvariation(spikes, freq, tselection, reg)

data = sort(spikes(:,2));
% Ntot = length(unique(data(:,1)));
difdata = diff(data);
dtj = mean(difdata);
dtj2 = mean(difdata.^2);
c = floor(tselection / (1000/mean(freq)));
J = length(difdata);
if c < 1 || J <= 2*(c+1)
    CV_all = NaN;
else
    CV_all = sqrt((c*(J-1))/(J*(J-c)))* (sqrt(dtj2 - dtj^2) / dtj);
end

CV_ind = zeros(reg,1);
for r = 1:reg
    data_ind = sort(spikes(spikes(:,3) == r, 2));
%     Na = length(unique(spikes(spikes(:,3) == r, 1)));
    difdata_ind = diff(data_ind);
    dtj_ind = mean(difdata_ind);
    dtj2_ind = mean(difdata_ind.^2);
    c = floor(tselection / (1000/freq(r)));
    J = length(difdata_ind);
    if c < 1 || J <= 2*(c+1)
        CV_ind(r) = NaN;
    else
        CV_ind(r) = sqrt((c*(J-1))/(J*(J-c)))* (sqrt(dtj2_ind - dtj_ind^2) / dtj_ind) ;
    end
end
end