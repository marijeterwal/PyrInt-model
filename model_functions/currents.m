function [Ii, Ie, Ii_dyn, Ie_dyn, ipulsevct, epulsevct] = currents(sett, para)

reg = sett.Nregions;
dt = sett.dt;
Ni = sett.Ni;
Ne = sett.Ne;
tnumb = sett.Ttot/dt +1;

if length(sett.loopIi) < reg || length(sett.loopIe) < reg || length(sett.iIu) < reg || length(sett.eIu) < reg;
    error('Currents settings are invalid: currents should be specified for all regions.')
end

%% Ie and Ii
Ii = zeros(sett.Ni*reg,1);
Ie = zeros(sett.Ne*reg,1);

for r = 1:reg
    s = RandStream('mcg16807','Seed',sett.seed+50+r); RandStream.setGlobalStream(s);
    Ii((r-1)*sett.Ni +1: r*sett.Ni,1) = normrnd(sett.iIu(r), sett.iIs, sett.Ni, 1);
    Ie((r-1)*sett.Ne +1: r*sett.Ne,1) = normrnd(sett.eIu(r), sett.eIs, sett.Ne, 1);
end

%% Time changing currents

Tseries = 0:dt:sett.Ttot;

% periodic input
if (sum(sett.Iiper) + sum(sett.Ieper) >= 1)
    tempT = Tseries - sett.Iper_phase;
    % change amp - constant offset
    Ii_per = single(sett.Iiper(1:reg)' * (sett.Iiper_offset + sett.Iper_amp * sin(2*pi*sett.Iper_freq*(tempT/1000) + sett.Iper_phase)));
    Ie_per = single(sett.Ieper(1:reg)' * (sett.Ieper_offset + sett.Iper_amp * sin(2*pi*sett.Iper_freq*(Tseries/1000) + sett.Iper_phase)));
%     % change offset - constant amp
%     Ii_per = single(sett.Iiper(1:reg)' * (sett.Iper_amp + 0.5 * sin(2*pi*sett.Iper_freq*(tempT/1000) + sett.Iper_phase)));
%     Ie_per = single(sett.Ieper(1:reg)' * (sett.Iper_amp + 0.5 * sin(2*pi*sett.Iper_freq*(Tseries/1000) + sett.Iper_phase)));
%     % change both
%     Ii_per = single(sett.Iiper(1:reg)' * (sett.Iper_amp + sett.Iper_amp * sin(2*pi*sett.Iper_freq*(tempT/1000) + sett.Iper_phase)));
%     Ie_per = single(sett.Ieper(1:reg)' * (sett.Iper_amp + sett.Iper_amp * sin(2*pi*sett.Iper_freq*(Tseries/1000) + sett.Iper_phase)));
    
%     Ii_per = single(sett.Iiper(1:reg)' * (sett.Iiper_offset + sett.Iper_amp/sett.Iper_ampfraction + sett.Iper_amp * sin(2*pi*sett.Iper_freq*(tempT/1000) + sett.Iper_phase))); % + Iper_amp/Iper_ampfraction
%     Ie_per = single(sett.Ieper(1:reg)' * (sett.Ieper_offset + sett.Iper_amp + sett.Iper_amp * sin(2*pi*sett.Iper_freq*(Tseries/1000) + sett.Iper_phase))); %  + Iper_amp
else
    Ii_per = single(zeros(reg, tnumb));
    Ie_per = single(zeros(reg, tnumb));
end

% ramp
if sum(sett.Iiramp) + sum(sett.Ieramp) >= 1;
    Ii_ramp = single(sett.Iiramp(1:reg)' * (sett.Iiramp_offset + ((sett.Iramp_max - sett.Iramp_min)* 1:tnumb)/ tnumb  / sett.Iramp_ampfraction));
    Ie_ramp = single(sett.Ieramp(1:reg)' * (sett.Ieramp_offset + ((sett.Iramp_max - sett.Iramp_min) * 1:tnumb)/tnumb));
else
    Ii_ramp = single(zeros(reg, tnumb));
    Ie_ramp = single(zeros(reg, tnumb));
end

% Noise
if sum(sett.Iinoise) + sum(sett.Ienoise) >= 1;
    s = RandStream('mcg16807','Seed',sett.noiseseed);
    RandStream.setGlobalStream(s);
    
    if sett.noisetype == 1 %White
        nsigw = randn(reg,tnumb);
        nsigw(reg,:) = sett.noisecor*nsigw(1,:) + (1-sett.noisecor)*nsigw(reg,:);
        Ii_noise = single(bsxfun(@times,sett.Iinoise(1:reg)',nsigw));
        Ie_noise = single(bsxfun(@times,sett.Ienoise(1:reg)',nsigw));
        clear nsigw
    elseif sett.noisetype == 2 %Brown
        alpha = 1-dt/sett.tau;
        afil = [1 -alpha];
        bfil = 1;
        nsigb = zeros(reg, tnumb);
        for r = 1:reg % because of 1D filter function
            randsig = sett.noiseamp * 0.5*randn(1,tnumb);
            nsigb(r,:) = filter(bfil,afil,randsig)*sqrt(1-alpha);
            nsigb(r,:) = sett.noiseamp*0.1 - mean(nsigb(r,:)) + nsigb(r,:);
        end
        nsigb(reg,:) = sett.noisecor*nsigb(1,:) + (1-sett.noisecor)*nsigb(reg,:);
        Ii_noise = single(bsxfun(@times,sett.Iinoise(1:reg)',nsigb));
        Ie_noise = single(bsxfun(@times,sett.Ienoise(1:reg)',nsigb));
        clear nsigb
    elseif sett.noisetype == 3 %Pink
        B = [0.049922035 -0.095993537 0.050612699 -0.004408786];
        A = [1 -2.494956002   2.017265875  -0.522189400];
        nsigp = zeros(reg, tnumb);
        for r = 1:reg
            randsig = sett.noiseamp * 2*randn(1,tnumb);
            nsigp(r,:) = filter(B,A,randsig);
        end
        nsigp(reg,:) = sett.noisecor*nsigp(1,:) + (1-sett.noisecor)*nsigp(reg,:);
        Ii_noise = single(bsxfun(@times,sett.Iinoise(1:reg)',nsigp)); 
        Ie_noise = single(bsxfun(@times,sett.Ienoise(1:reg)',nsigp));
        clear nsigp
    end
else
    Ii_noise = single(zeros(reg, tnumb));
    Ie_noise = single(zeros(reg, tnumb));
end

% Step currents
if sum(sett.Iistep) + sum(sett.Iestep) > 0
               
    % check definitions of amplitude
    l_amp = arrayfun(@(x) length(sett.Iestep_amp{x}), 1:length(sett.Iestep_amp));
    assert(sum(diff(l_amp)==0)>0, 'Error: amplitudes of step currents are not of equal lengths')
    l_amp = arrayfun(@(x) length(sett.Iistep_amp{x}), 1:length(sett.Iistep_amp));
    assert(sum(diff(l_amp)==0)>0, 'Error: amplitudes of step currents are not of equal lengths')
    nsteps = max(l_amp);
    
    % check length of run
    l_time = nsteps*sett.Istep_dur/dt + sett.Tselection/dt + 1;
    assert(l_time == tnumb, 'Length of run does not agree with step current definitions');
    
    Ii_step = zeros(reg, tnumb, 'single');
    Ie_step = zeros(reg, tnumb, 'single');  
    
    
    for r = 1:reg
        Ii_step(r,1:sett.Tselection/dt) = sett.Iistep(r)*sett.Iistep_amp{r}(1);
        Ie_step(r,1:sett.Tselection/dt) = sett.Iestep(r)*sett.Iestep_amp{r}(1);
        for ns = 1:nsteps
            Ii_step(r,sett.Tselection/dt + (ns-1)*sett.Istep_dur/dt+1:sett.Tselection/dt + (ns)*sett.Istep_dur/dt) = ...
                sett.Iistep(r)*sett.Iistep_amp{r}(ns);
            Ie_step(r,sett.Tselection/dt + (ns-1)*sett.Istep_dur/dt+1:sett.Tselection/dt + (ns)*sett.Istep_dur/dt) = ...
                sett.Iestep(r)*sett.Iestep_amp{r}(ns);
        end
    end
    
else
    Ii_step = zeros(reg, tnumb, 'single');
    Ie_step = zeros(reg, tnumb, 'single'); 
end

if sett.Istep_repNoise
    Ii_noise = Ii_noise(:,1:(sett.Tselection+sett.Istep_dur)/dt+1);
    Ie_noise = Ie_noise(:,1:(sett.Tselection+sett.Istep_dur)/dt+1);
    for ns = 2:nsteps
        Ii_noise = cat(2,Ii_noise, Ii_noise(:,sett.Tselection/dt+1:(sett.Tselection+sett.Istep_dur)/dt));
        Ie_noise = cat(2,Ie_noise, Ie_noise(:,sett.Tselection/dt+1:(sett.Tselection+sett.Istep_dur)/dt));
    end    
end


%% determine which neurons receive the currents
Ii_dyn_vect = zeros(Ni*reg, 1); Ii_dyn_vect(rand(Ni*reg,1) <= sett.Idyn_frac) = 1;
Ie_dyn_vect = zeros(Ne*reg, 1); Ie_dyn_vect(rand(Ne*reg,1) <= sett.Idyn_frac) = 1;

for r = 1:reg
    if sett.subregion_noise == 0 && sett.subregion >= r
        Ii_dyn_vect(r*Ni-Nisub+1:r*Ni) = 0;
        Ie_dyn_vect((r-1)*Ne+1:(r-1)*Ne+Nesub) = 0;
    elseif sett.subregion_noise == 2 && sett.subregion >= r
        Ii_dyn_vect((r-1)*Ni+1:r*Ni-Nisub) = 0;
        Ie_dyn_vect((r-1)*Ne+Nesub+1:r*Ne) = 0;
    end

    Ii_dyn((r-1)*Ni+1:r*Ni,:) = Ii_dyn_vect((r-1)*Ni+1:r*Ni) * (Ii_per(r,:) + Ii_ramp(r,:) + Ii_noise(r,:) + Ii_step(r,:));
    Ie_dyn((r-1)*Ne+1:r*Ne,:) = Ie_dyn_vect((r-1)*Ne+1:r*Ne) * (Ie_per(r,:) + Ie_ramp(r,:) + Ie_noise(r,:) + Ie_step(r,:));
end

clear Ii_per Ii_ramp Ii_noise Ie_per Ie_ramp Ie_noise Ii_step Ie_step

%% Pulses

tnumbpulse = 30/dt+1;
ipulsevct = zeros(Ni*reg,tnumbpulse);
epulsevct = zeros(Ne*reg,tnumbpulse);

if sum(sett.Iipulse) + sum(sett.Iepulse) >= 1;
    
    if sett.pulse_amptype == 1;
        
        g_syn_ipulse = sett.g_syn_ipulse;
        g_syn_epulse = sett.g_syn_epulse;
        ealpha_syn = para.ealpha_syn;
        ebeta_syn = para.ebeta_syn;
        
        % amp = percentage of neurons stimulated
        p = sett.Ipulse_amp/100;
        
        s = RandStream('mcg16807','Seed',sett.seed+20); RandStream.setGlobalStream(s);
        pulsevectori = zeros(Ni*reg, 1);
        pulsevectori(rand(Ni*reg, 1) < p) = 1;
        
        s = RandStream('mcg16807','Seed',sett.seed+25); RandStream.setGlobalStream(s);
        pulsevectore = zeros(Ne*reg, 1);
        pulsevectore(rand(Ne*reg, 1) < p) = 1;
                
        t1 = 0:dt:sett.Ipulse_length;
        t2 = (sett.Ipulse_length+dt):dt:30;
        
        if ~isempty(find(sett.Iipulse,reg)) || ~isempty(find(sett.Iepulse,reg))
            for r = 1:reg
                if sett.Iipulse(r) == 1;
                    ipulsevct((r-1)*Ni + 1: r*Ni, 1:length(t1)) = pulsevectori((r-1)*Ni + 1: r*Ni) * g_syn_ipulse * (-(ealpha_syn/(ealpha_syn + ebeta_syn))*exp(-(ealpha_syn + ebeta_syn)*t1) + (ealpha_syn/(ealpha_syn + ebeta_syn))) ;
                    ipulsevct((r-1)*Ni + 1: r*Ni,length(t1)+1:end) = pulsevectori((r-1)*Ni + 1: r*Ni) * g_syn_ipulse * (ealpha_syn/(ealpha_syn + ebeta_syn))* (exp(ebeta_syn)-exp(-ealpha_syn))*exp(- ebeta_syn*t2);
                end
                if sett.Iepulse(r) == 1
                    epulsevct((r-1)*Ne + 1: r*Ne,1:length(t1)) = pulsevectore((r-1)*Ne + 1: r*Ne) * g_syn_epulse * (-(ealpha_syn/(ealpha_syn +ebeta_syn))*exp(-(ealpha_syn + ebeta_syn)*t1) + (ealpha_syn/(ealpha_syn + ebeta_syn)));
                    epulsevct((r-1)*Ne + 1: r*Ne,length(t1)+1:end) = pulsevectore((r-1)*Ne + 1: r*Ne) * g_syn_epulse * (ealpha_syn/(ealpha_syn + ebeta_syn))* (exp(ebeta_syn)-exp(-ealpha_syn))*exp(- ebeta_syn*t2);
                end
            end
        end
        
    elseif sett.pulse_amptype == 2;
        % pulses - amp = strength of gsyn_pulse
        p = sett.Ipulse_amp/100;
        pulsevectori = ones(Ni*reg, 1);
        pulsevectore = ones(Ne*reg, 1);
        
        g_syn_ipulse = p * sett.g_syn_ipulse;
        g_syn_epulse = p * sett.g_syn_epulse;
        ealpha_syn = para.ealpha_syn;
        ebeta_syn = para.ebeta_syn;
        
        t1 = 0:dt:sett.Ipulse_length;
        t2 = (sett.Ipulse_length+dt):dt:30;
        
        if ~isempty(find(ipulse,reg)) || ~isempty(find(epulse,reg))
            for r = 1:reg
                if ipulse(r) == 1;
                    ipulsevct((r-1)*Ni + 1: r*Ni, 1:length(t1)) = pulsevectori((r-1)*Ni + 1: r*Ni) * g_syn_ipulse * (-(ealpha_syn/(ealpha_syn + ebeta_syn))*exp(-(ealpha_syn + ebeta_syn)*t1) + (ealpha_syn/(ealpha_syn + ebeta_syn))) ;
                    ipulsevct((r-1)*Ni + 1: r*Ni,length(t1)+1:end) = pulsevectori((r-1)*Ni + 1: r*Ni) * g_syn_ipulse * (ealpha_syn/(ealpha_syn + ebeta_syn))* (exp(ebeta_syn)-exp(-ealpha_syn))*exp(- ebeta_syn*t2);
                end
                if epulse(r) == 1
                    epulsevct((r-1)*Ne + 1: r*Ne,1:length(t1)) = pulsevectore((r-1)*Ne + 1: r*Ne) * g_syn_epulse * (-(ealpha_syn/(ealpha_syn +ebeta_syn))*exp(-(ealpha_syn + ebeta_syn)*t1) + (ealpha_syn/(ealpha_syn + ebeta_syn)));
                    epulsevct((r-1)*Ne + 1: r*Ne,length(t1)+1:end) = pulsevectore((r-1)*Ne + 1: r*Ne) * g_syn_epulse * (ealpha_syn/(ealpha_syn + ebeta_syn))* (exp(ebeta_syn)-exp(-ealpha_syn))*exp(- ebeta_syn*t2);
                end
            end
        end
    end
end
