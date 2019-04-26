% Standard analyses of measures between areas   

numbercomb = nchoosek(reg,2);

comb = flipud(combnk(1:reg,2));

%% synchrony, delay and phase difference
phaseosc1 = phaseosc(1,Tselection/dthist:end-100,:);
phaseosci1 = phaseosci(1,Tselection/dthist:end-100,:);

for nc = 1:numbercomb
    reg1 = comb(nc,1);
    reg2 = comb(nc,2);
    
%     [cor12, ~] = crosscorrelation(hist_e(:,Tselection/dthist:end,reg1), hist_e(:,Tselection/dthist:end,reg2), []);
    cor12 = xcorr(hist_e(:,Tselection/dthist:end,reg1), hist_e(:,Tselection/dthist:end,reg2), 'coeff');
    dum = findpeaks(cor12, 'sortstr', 'descend', 'npeaks', 1);
    if isempty(dum)
        peakcor(l1,l2,nc) = NaN;
    else
        peakcor(l1,l2,nc) = dum;
    end
%     [peakcor(l1,l2,nc), ~] = determinecorpeak(cor12, dthist);
    
    phaseshiftseries =  wrapTo2Pi((unwrap(phaseosc1(:,:,reg1)*2*pi) - unwrap(phaseosc1(:,:,reg2)*2*pi)))/ (2*pi);
    phasediff(l1,l2,nc) = mean(phaseshiftseries);
    phaseshiftseriesi =  wrapTo2Pi((unwrap(phaseosci1(:,:,reg1)*2*pi) - unwrap(phaseosci1(:,:,reg2)*2*pi)))/ (2*pi);
    phasediffi(l1,l2,nc) = mean(phaseshiftseriesi);
    phasediff_std(l1,l2,nc) = std(phaseshiftseries);
    
    delaye(l1,l2,nc) = drivedelay(cor12, dthist); %(round(length(corei)/2):end)
    delay(l1,l2,nc) = delaye(l1,l2,nc);
    if delay(l1,l2,nc) < -1000/avgmorfreq(l1,l2,2); delay(l1,l2,nc) = NaN;
    elseif delay(l1,l2,nc) < 0; delay(l1,l2,nc) = 1000./avgmorfreq(l1,l2,2) + delay(l1,l2,nc); 
    end

    % coherence
    [cohe, phase,~,~,~,f,~,~,~] = coherencyc(hist_e(1,Tselection/dthist:length(hist_e(1,:,reg1)),reg1),hist_e(1,Tselection/dthist:length(hist_e(1,:,reg1)),reg2), params);
%     [cohi, ~,~,psdi1,psdi2,~,~,~,~] = coherencyc(hist_i(1,Tselection/dthist:length(hist_i(1,:,reg1)),reg1),hist_i(1,Tselection/dthist:length(hist_i(1,:,reg1)),reg2), params);
%     psd_e = [psde1,psde2];
%     psd_i = [psdi1,psdi2];
%     fpsd = f;
    fcoh = f;
    [~, fId] = min(abs(fcoh - avgmorfreq(l1,l2,reg1)));
    coh(l1,l2,nc) = cohe(fId);
    
    [cohit, phasi,~,~,~,f,~,~,~] = coherencyc(hist_i(1,Tselection/dthist:length(hist_i(1,:,reg1)),reg1),hist_i(1,Tselection/dthist:length(hist_i(1,:,reg1)),reg2), params);
    fcoh = f;
    [~, fId] = min(abs(fcoh - avgmorfreq(l1,l2,reg1)));
    cohi(l1,l2,nc) = cohit(fId);
    
%     [Cxy, fCxy] = mscohere(hist_e(1,Tselection/dthist:length(hist_e(1,:,reg1)),reg1),...
%         hist_e(1,Tselection/dthist:length(hist_e(1,:,reg2)),reg2), [], [], [], params.Fs);
%     [~, fId] = min(abs(fCxy - avgmorfreq(l1,l2,reg1)));
%     Ccoh(l1,l2,nc) = Cxy(fId);
%     
%     psd_e(:,1) = pmtm(hist_e(1,Tselection/dthist:length(hist_e(1,:,reg1)),reg1), 1.25);
%     psd_e(:,2) = pmtm(hist_e(1,Tselection/dthist:length(hist_e(1,:,reg2)),reg2), 1.25);

    %% phase locking - spikes of area 2 with lfp of area 1
    
    if analyzePPC
        combnU = [nchoosek(1:reg,2); fliplr(nchoosek(1:reg,2))];
        %
        Phaselocki_12(l1,l2,:) = phaselocking(sett, spikesselection_i(spikesselection_i(:,3) == 2,:), phaseosc(1,:,1), [2,1]);
        Phaselocke_12(l1,l2,:) = phaselocking(sett, spikesselection_e(spikesselection_e(:,3) == 2,:), phaseosc(1,:,1), [2,1]);
    elseif analyzePPC == 0 && loadanalysis == 0 && l1 == l1steps && l2 == l2steps
        load([sett.savelocdata,'analysis2', '.mat'], 'Phaselocki_12')
        load([sett.savelocdata,'analysis2', '.mat'], 'Phaselocke_12')
    end
    
    %% old
%     Phaselocki_12(l1,l2,:) = phaselocking(sett, spikesselection_i(spikesselection_i(:,3) == 2,:), phaseosc(1,:,1), Ni);
%     Phaselocke_12(l1,l2,:) = phaselocking(sett, spikesselection_e(spikesselection_e(:,3) == 2,:), phaseosc(1,:,1), Ne);
    
%     Pl_ni = zeros(Ni,1);
%     for nn = 1:Ni
%         spikesNr = round(spikesselection_i(spikesselection_i(:,1) == nn & spikesselection_i(:,3) == reg2,2)/dthist);
%         phaseNr = phaseosc(1,spikesNr,reg1);
%         %    phaseNr = mod(spikesNr, mean(diff(spikesNr))) / mean(diff(spikesNr));
%         Pl_ni(nn) = abs(sum(exp(1i*phaseNr*2*pi))/Ni);
%     end
%     
%     Pl_ne = zeros(Ne,1);
%     for nn = 1:Ne
%         spikesNr = round(spikesselection_e(spikesselection_e(:,1) == nn & spikesselection_e(:,3) == reg2,2)/dthist);
%         phaseNr = phaseosc(1,spikesNr, reg1);
%         Pl_ne(nn) = abs(sum(exp(1i*phaseNr*2*pi))/Ne);
%     end
%     Phaselocke_12(l1,l2,nc) = mean(Pl_ne);


end

        %

