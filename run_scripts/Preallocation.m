% preallocate memory for the variables used to save analyses

if reg>1
    ncomb = nchoosek(reg,2);
else 
    ncomb = 1;
end

%% standard analysis for each region

freqi = zeros(l1steps, l2steps, reg);
freqe = zeros(l1steps, l2steps, reg);
avgmorfreq = zeros(l1steps, l2steps, reg);
avgmorfreqi = zeros(l1steps, l2steps, reg);
CVi = zeros(l1steps, l2steps, reg);
CVe = zeros(l1steps, l2steps, reg);
CVall = zeros(l1steps, l2steps, reg);
ri = zeros(l1steps, l2steps, reg);
re = zeros(l1steps, l2steps, reg);
delayei = zeros(l1steps, l2steps, reg);
spikesnumbers_e1 = zeros(l1steps, l2steps, reg);
spikesnumbers_e2 = zeros(l1steps, l2steps, reg);
Phaselocki = zeros(l1steps, l2steps, reg);
Phaselocke = zeros(l1steps, l2steps, reg);

varNames_standard1  = {'freqi', 'freqe', 'avgmorfreq','avgmorfreqi', 'CVi', 'CVe', 'CVall', ...
    'ri', 're', 'delayei', 'Phaselocki', 'Phaselocke'};

%% standard analysis for two or more regions

if reg > 1
    peakcor = zeros(l1steps, l2steps, ncomb);
    coh = zeros(l1steps, l2steps, ncomb);
    cohi = zeros(l1steps, l2steps, ncomb);
    psd_e = zeros(l1steps, l2steps, ncomb);
    delay = zeros(l1steps, l2steps, ncomb);
    delaye = zeros(l1steps, l2steps, ncomb);
    phasediff = zeros(l1steps, l2steps, ncomb);
    phasediffi = zeros(l1steps, l2steps, ncomb);
    phasediff_std = zeros(l1steps, l2steps, ncomb);
    Phaselocki_12 = zeros(l1steps, l2steps, ncomb);
    Phaselocke_12 = zeros(l1steps, l2steps, ncomb);
    
    varNames_standard2  = {'peakcor', 'coh', 'cohi', 'delay', 'phasediff', 'phasediffi', 'phasediff_std',...
    'Phaselocki_12', 'Phaselocke_12'};
end


%% loop variable specific analyses

if loopId(5) > 0
    
    phaseOfPulse = zeros(l1steps, l2steps, reg);
    
    local_phaseshift = zeros(l1steps, l2steps, reg);
    local_phaseshift_full = zeros(l1steps, l2steps, reg);
    local_phaseshift_std = zeros(l1steps, l2steps, reg);
    
    reg_phaseshift = zeros(l1steps, l2steps, ncomb);
    reg_phaseshift_full = zeros(l1steps, l2steps, ncomb);
    reg_phaseshift_std = zeros(l1steps, l2steps, ncomb);
    
%     if analyzeInfo == true
        h_base = zeros(l1steps, l2steps, reg);
        h_fper = zeros(l1steps, l2steps, reg);
        h_spikes = zeros(l1steps, l2steps, reg);
        h_freq = zeros(l1steps, l2steps, reg);
        h_power = zeros(l1steps, l2steps, reg);
        h_peak = zeros(l1steps, l2steps, reg);
        h_peakheight = zeros(l1steps, l2steps, reg);
        h_peaktime = zeros(l1steps, l2steps, reg);
        h_sync = zeros(l1steps, l2steps, reg);
        
        deltaTimes_First_Input = zeros(l1steps, l2steps, reg);
        deltaHeights_First_Input = zeros(l1steps, l2steps, reg);
        deltaTimes = zeros(l1steps, l2steps, reg, npeaks);
        deltaHeights = zeros(l1steps, l2steps, reg, npeaks);
        
        varNames_pulses = {'phaseOfPulse', ...
            'h_base', 'h_fper', 'h_spikes', 'h_freq', 'h_power', 'h_peak', ...
            'h_peaktime', 'h_peakheight', 'h_sync',...
            'deltaTimes_First_Input', 'deltaHeights_First_Input', ...
            'deltaTimes', 'deltaHeights'};
%     else
%         varNames_pulses = {'phaseOfPulse'};
%     end
end


%% selected analyses

if analyzeInfo == true
    nMI_freq    = zeros(l1steps, l2steps, 2*reg-1);
    nMI_peak    = zeros(l1steps, l2steps, 2*reg-1);
    nMI_power   = zeros(l1steps, l2steps, 2*reg-1);
    nMI_base    = zeros(l1steps, l2steps, 2*reg-1);
    nMI_fper    = zeros(l1steps, l2steps, 2*reg-1);
    
    rho_freq    = zeros(l1steps, l2steps, 2*reg-1);
    p_freq      = zeros(l1steps, l2steps, 2*reg-1);
    rho_peak    = zeros(l1steps, l2steps, 2*reg-1);
    p_peak      = zeros(l1steps, l2steps, 2*reg-1);
    rho_base    = zeros(l1steps, l2steps, 2*reg-1);
    p_base      = zeros(l1steps, l2steps, 2*reg-1);
    rho_power   = zeros(l1steps, l2steps, 2*reg-1);
    p_power     = zeros(l1steps, l2steps, 2*reg-1);
    rho_fper   = zeros(l1steps, l2steps, 2*reg-1);
    p_fper     = zeros(l1steps, l2steps, 2*reg-1);
    
    MI_freq     = zeros(l1steps, l2steps);
    MI_base     = zeros(l1steps, l2steps);
    MI_power    = zeros(l1steps, l2steps);
    MI_fper    = zeros(l1steps, l2steps);
    
    rho_freq_12 = zeros(l1steps, l2steps, ncomb);
    rho_base_12 = zeros(l1steps, l2steps, ncomb);
    rho_power_12 = zeros(l1steps, l2steps, ncomb);
    rho_peak_12 = zeros(l1steps, l2steps, ncomb);
    rho_fper_12 = zeros(l1steps, l2steps, ncomb);
    
    varNames_info = {'rho_freq', 'rho_peak', 'rho_base', 'rho_power','rho_fper', ...
        'p_freq', 'p_peak', 'p_base', 'p_power', ...
        'rho_freq_12', 'rho_base_12', 'rho_power_12', 'rho_peak_12', 'rho_fper_12',...
            'MI_freq', 'MI_power', 'MI_base', 'MI_fper',...
            'nMI_freq', 'nMI_peak', 'nMI_power', 'nMI_base','nMI_fper'};
end
