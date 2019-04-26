% Standard analyses for each area

 % not valid for subregions

%% histograms

[hist_i, histtime_i, smoothhist_i] = histograms(spikes_i, Ni, Ttot, dt, dthist, cutoff_freq_low, cutoff_freq_high, reg);
[hist_e, histtime_e, smoothhist_e] = histograms(spikes_e, Ne, Ttot, dt, dthist, cutoff_freq_low, cutoff_freq_high, reg);

%% analysis on selection of the STH
spikesselection_i = spikes_i((spikes_i(:,2,:) >= Tselection), :);
spikesselection_e = spikes_e((spikes_e(:,2,:) >= Tselection), :);
selhisti = smoothhist_i(1,Tselection/dthist:end,:);
selhiste = smoothhist_e(1,Tselection/dthist:end,:);

spikesnumbers_e1(l1,l2,1) = length(find(spikesselection_e(:,3)==1));
spikesnumbers_e2(l1,l2,1) = length(find(spikesselection_e(:,3)==2));

% firing rate
ri(l1,l2,:) = firingrate(spikesselection_i,Ni,reg);
re(l1,l2,:) = firingrate(spikesselection_e,Ne,reg);

% oscillation frequency based on auto-correlation
[corii, ~] = autocorrelation(selhisti, corselection/dthist, reg);
[coree, ~] = autocorrelation(selhiste, corselection/dthist, reg);
freqi(l1,l2,:) = osc_freq(corii, dthist, reg);
freqe(l1,l2,:) = osc_freq(coree, dthist, reg);


% [freqosc1, phaseosc1, power1, ~] = freq_phase(hist_e(1,1:3200/dthist,:), dthist, sett.corselection, reg);
% [freqosc2, phaseosc2, power2, ~] = freq_phase(hist_e(1,3200/dthist:6200/dthist,:), dthist, sett.corselection, reg);
% [freqosc3, phaseosc3, power3, ~] = freq_phase(hist_e(1,6200/dthist:end,:), dthist, sett.corselection, reg);
% freqosc = cat(2,freqosc1,freqosc2,freqosc3);
% phaseosc = cat(2,phaseosc1,phaseosc2(:,2:end,:),phaseosc3(:,2:end,:));
% power = cat(2,power1,power2(:,2:end,:),power3(:,2:end,:));

[freqosc, phaseosc, power, powermat] = freq_phase(hist_e, dthist, sett.corselection, reg, []);
[freqosci, phaseosci, poweri, powermati] = freq_phase(hist_i, dthist, sett.corselection, reg, []);
avgmorfreq(l1,l2,:) = nanmean(freqosc(1,Tselection/dthist:end-100,:), 2);
avgmorfreqi(l1,l2,:) = nanmean(freqosci(1,Tselection/dthist:end-100,:), 2);

% coefficient of variation
[CVi(l1,l2,:), ~] = coefficientofvariation(spikesselection_i, freqi(l1,l2,:), Ttot - Tselection, reg); % CV_ind is a vector with the CVp of the individual regions
[CVe(l1,l2,:), ~] = coefficientofvariation(spikesselection_e, freqe(l1,l2,:), Ttot - Tselection, reg);
[CVall(l1,l2,:), ~] = coefficientofvariation([spikesselection_e; spikesselection_i], freqe(l1,l2,:), Ttot - Tselection, reg);

for r = 1:reg
    corei = xcorr(hist_e(:,Tselection/dthist:end,r), hist_i(:,Tselection/dthist:end,r), 'coeff');
    delayei(l1,l2,r) = drivedelay(corei, dthist); %(round(length(corei)/2):end)
    if delayei(l1,l2,r) < -1000/avgmorfreq(l1,l2,r); delayei(l1,l2,r) = NaN;
    elseif delayei(l1,l2,r) < 0; delayei(l1,l2,r) = 1000./avgmorfreq(l1,l2,r) + delayei(l1,l2,r); 
    end
end
         
%% phase-locking

if analyzePPC
    Phaselocki(l1,l2,:) = phaselocking(sett, spikesselection_i, phaseosc, []);
    Phaselocke(l1,l2,:) = phaselocking(sett, spikesselection_e, phaseosc, []);
elseif analyzePPC == 0 && loadanalysis == 0 && l1 == l1steps && l2 == l2steps
    load([sett.savelocdata,'analysis', '.mat'], 'Phaselocki') % safety: this is an expensive analysis. Don't overwrite with zeros
    load([sett.savelocdata,'analysis', '.mat'], 'Phaselocke')
end

if sum(sett.Ieper)>0 && reg == 1
    [~, ~, Ii_dyn, Ie_dyn, ~, ~] = currents(sett, para);
    phaselockOsce(l1,l2,:) = phaselocking(sett, spikes_e, Ie_dyn, []);
    phaselockOsci(l1,l2,:) = phaselocking(sett, spikes_i, Ie_dyn, []);
end
