function [sett, para] = settings(savesettings, runnumber)
sett = struct;
para = struct;

sett.runnumber = runnumber;

%% time steps, population size and initial conditions
sett.dt = 0.05; %ms, time step
sett.Ttot = 1200; % ms, duration of simulation

sett.Ne = 400; 
sett.Ni = 100;

sett.Nregions = 2;

sett.fixed_initcond = 2; % 1 == identical for all I (-65) and all E (-75 mV); 2 == rand around -50 mV;
sett.seed = 2; % for restarting random seed generator

%% analysis
sett.Tselection = 200; % ms
sett.corselection = 100; % ms, maximum delay for autocorrelation

% histograms
sett.dthist = 0.5;% ms
sett.cutoff_freq_low = 10;
sett.cutoff_freq_high = 200;

%% Connections
sett.connectiontype = 3; %1 = all-to-all (Msyn not required),
                        %2 = no connections
                        %3 = unique Msyn values for different connection types, 
sett.networktype = 3;    %1 = all-to-all (Msyn not required),
                        %2 = no connections
                        
                        %3 = unique Msyn values for different connection types,
% connection probabilities between cell types within region
sett.Msynii = 0.40; % of i cells send connections to each i cell, i.e convergence
sett.Msynie = 0.60; % of i cells send connections to each e cell
sett.Msynee = 0.10; 
sett.Msynei = 0.20;

% connection probabilities between cell types between regions
sett.Msynee_rtor = [[0.00, 0.05]; % this should be a reg*reg matrix; diagonals should be 0; 
                    [0.05, 0.00]]; % row = from; column = to
sett.Msynei_rtor = [[0.00, 0.10]; % this should be a reg*reg matrix; diagonals should be 0;
                    [0.10, 0.00]]; % row = from; column = to
sett.Msynii_rtor = [[0.00, 0.10]; % this should be a reg*reg matrix; diagonals should be 0;
                    [0.10, 0.00]]; % row = from; column = to
sett.Msynie_rtor = [[0.00, 0.10]; % this should be a reg*reg matrix; diagonals should be 0;
                    [0.10, 0.00]]; % row = from; column = to           
%% driving currents

% loop over driving currents
sett.loopIi = [5,0]; % 0 is no loop, 1 = loop 5 = specialloop (see later);
sett.loopIe = [5,0];
sett.Istep = 0.1;

sett.iIu = [1.0,1.0]; %uA/cm^2, average amplitude of applied current; if reg > 1 and length(iIu) == 1, this value is used for all regions, if reg ~= length(iIu) an error occurs
sett.iIu_max = [4,5];

sett.eIu = [1.5,1.5]; %uA/cm^2
sett.eIu_max = [2,1.7173];

sett.iIspecial{1} = linspace(0.4, 2.0, 40);
% Sync at low gamma: linspace(1.1, 0.86, 30);
% Freq at low sync (0.06): linspace(0.5, 2.6, 50);
% Freq at med sync (0.085): linspace(0.4, 2.0, 40);
% Freq at high sync (0.125): linspace(0.6,1.2,30);
sett.eIspecial{1} = 0.4652*sett.iIspecial{1}.^4 - 1.9860*sett.iIspecial{1}.^3 + 3.2879*sett.iIspecial{1}.^2 - 1.0623*sett.iIspecial{1} + 0.7546;
% Sync at low gamma: -0.5*sett.iIspecial + 2;
% Freq at low sync: 0.0983*sett.iIspecial.^3 + -0.1679*sett.iIspecial.^2 + 1.3415*sett.iIspecial + 0.0961;
% Freq at med sync: 0.4652*sett.iIspecial.^4 - 1.9860*sett.iIspecial.^3 + 3.2879*sett.iIspecial.^2 - 1.0623*sett.iIspecial + 0.7546;
% Freq at high sync: 2.0700*sett.iIspecial.^2 + -2.1058*sett.iIspecial + 1.6995;

sett.lateonset = 0; % late onset of drive to reg 2
sett.lateonset_start = 20; %16
sett.lateonset_end = 20;%184; %ms
sett.lateonset_steps = 1; % 1 %ms
                    
% heterogeneity of currents
sett.iIs = 0.1; % uA/cm^2 stdev of applied current; 0.02 low; 0.1 high
sett.eIs = 0.1; %uA/cm^2; 0.1 low; 0.3 high
% noise
sett.noise = 'low'; % 'very low, 'low', 'high' or 'no


%% applied currents

% fraction of cells affected:
sett.Idyn_frac = 1; % not for pulses! 

% noisy input per region
sett.Iinoise = [0,0];%0 = no noise; 1 = noise;
sett.Ienoise = [0,0];
sett.Inoise_loop = 0; % change noise amp

sett.noisetype = 2; % 1 = white noise; 2 = brownian noise; 3 = Pink;
sett.tau = 200; %in ms, for brownian noise
sett.subregion_noise = 0; %0: just into region; 1 = into region and subregion; 2 = just into subregion
sett.noiseamp = 0.5;
sett.noiseamp_max = 1.5;
sett.noiseamp_step = 1;%
sett.noisecor = 0; % creates a correlation between area 1 and area reg
sett.noiseseed = sett.seed+ 19; 

% pulses
sett.Iipulse = [0,0]; % 0 = no applied current; 1 = pulse;
sett.Iepulse = [0,0];
sett.pulse_amptype = 1; %1 = #neurons; 2 = gsyn;

sett.Iapploop = 0; % 0 = no loop; 1 = loop over start times; 2 = loop over amplitudes; 3 = loop over both amp and time
sett.timestep = 2; % ms
sett.ampstep = 2.5; % 

sett.Ipulse_start = 340; % ms
sett.Ipulse_end = 380;
sett.Ipulse_length = 1; % ms
sett.Ipulse_amp = 75; 
sett.Ipulse_amp_max = 100;

% sinusoid
sett.Iiper = [0,0]; % 0 = no applied current; 1 = oscillatory input to region;
sett.Ieper = [0,0];

sett.Iperloop = 0; % 0 = no loop; 1 = loop over periods; 2 loop over amplitude ; 3 is loop over amp & periods
sett.perampstep = -0.1;
sett.perstep = 5; %hz

% sett.Iper_period = 10; %ms
% sett.Iper_period_max = 50; %ms
sett.Iper_freq = 25; %freq
sett.Iper_freq_max = 25; %freq
sett.Iper_amp = 0; %uA
sett.Iper_amp_max = -2.5; %uA
sett.Iiper_offset = 0;
sett.Ieper_offset = -1;%-sett.Iper_amp/2;
sett.Iper_ampfraction = 1;
sett.Iper_phase = 0; %ms

% ramp
sett.Iiramp = [0,0];
sett.Ieramp = [0,0];

sett.Iramp_min = 0;
sett.Iramp_max = 1;
sett.Ieramp_offset = 0;
sett.Iiramp_offset = 0;
sett.Iramp_ampfraction = 2.5;

% stepwise changing currents
sett.Iistep = [0,0];
sett.Iestep = [0,0];

sett.Istep_dur = 3000; % step duration - ms
sett.Iistep_amp = {[0,0,0],[0,0,0],[0,0,0]};%
sett.Iestep_amp = {[0,0,0],[0,0,0],[0,0,0]};%
sett.Istep_repNoise = 0;%1; % 0 or 1

%% connections 
%values are total gsyns; divide by Njj to obtain unitary gsyn. Values are
%given after % for connectiontype == 3
sett.g_syn_ii = 0.0120 * sett.Msynii*sett.Ni; %1.1; %mS/cm2
sett.g_syn_ie = 0.0050 * sett.Msynie*sett.Ni; %0.24; 
sett.g_syn_ee = 0.0012 * sett.Msynee*sett.Ne;  %0.01; 
sett.g_syn_ei = 0.0010 * sett.Msynei*sett.Ne; %0.34; 

sett.g_syn_ee_r = [[0.00, 0.00];
                   [0.02, 0.00]];
sett.g_syn_ei_r = [[0.00, 0.00];
                   [0.02, 0.00]];
sett.g_syn_ii_r = [[0.00, 0.00];
                   [0.00, 0.00]];
sett.g_syn_ie_r = [[0.00, 0.00];
                   [0.00, 0.00]];

sett.g_syn_loop = [[0, 0]; % 1 is loop for E projections; -1 for I projections; 2 is for both
                   [1, 0]];

sett.g_syn_r_max = 0.30;
sett.g_syn_r_step = 0.02;
sett.g_syn_ifactor = 1;

%for pulse, gsyn values are unitary
sett.g_syn_epulse = 0.10;
sett.g_syn_ipulse = 0.10;

sett.delay = 1; %ms
sett.delay_reg_e = 5; %ms
sett.delay_reg_i = 5; %ms requirement: delay_reg_i >= delay_reg_e

%% Testing communication
sett.Icom = 0; %0 is no com. signal; 1 = into (sub of) region1;

sett.Iicom = [0,0]; % only works if sett.subregion_type = 1 or 2
sett.Iecom = [0,0];

sett.com_numneuron = 50; % percentage of cells
sett.comtype = 1; %1 = spikes, 2 = periodic, 3 = periodically modulated pulse train

sett.com_rate = 2; %Hz
sett.com_amp = 0;%only relevant if comtype == 2

% Subregions
sett.subregion = 0; %0 = no subregion; 1 = subregion in region 1; 2 = subregions in region 1 and 2

%% noise - continued
if strcmp(sett.noise,'low');
    sett.iIlambda = 0.02; %mV^2/ms; low = 0.02; high = 0.2
    sett.eIlambda = 0.06; %mV^2/ms; low = 0.1; high = 0.6
elseif strcmp(sett.noise,'very low');
    sett.iIlambda = 0.001; %mV^2/ms; low = 0.02; high = 0.2
    sett.eIlambda = 0.001; %mV^2/ms; low = 0.1; high = 0.6
elseif strcmp(sett.noise,'high');
    sett.iIlambda = 0.2; %mV^2/ms; low = 0.02; high = 0.2
    sett.eIlambda = 0.6; %mV^2/ms; low = 0.1; high = 0.6
else
    sett.iIlambda = 0; %mV^2/ms; low = 0.02; high = 0.2
    sett.eIlambda = 0; %mV^2/ms; low = 0.1; high = 0.6
end

%% parameters of HH cells (loaded in calcEInetwork functions)
para.iphi = 5;

%Reverse potentials
para.iE_l = -65; %mV
para.iE_na = 55;
para.iE_k = -90;

para.eE_l = -70; %mV
para.eE_na = 55;
para.eE_k = -90;

para.E_ampa = 0; % mV
para.E_gaba = -75;

%Conductances
para.ig_l = 0.1; %mS/cm^2
para.ig_na = 35;
para.ig_k = 9;
para.eg_l = 0.02; %mS/cm^2
para.eg_na = 24;
para.eg_nap = 0.07;
para.eg_kdr = 3;
para.eg_ka = 1.4;
para.eg_kslow = 1;

%Synapses
para.itheta_syn = 0; %mV
para.etheta_syn = -20; %mV
para.ialpha_syn = 10; % ms^-1; Bartos: 6.25 ; Wang: 12
para.ibeta_syn = 0.2; %ms^-1; decay time of synaptic conductance = 1/beta_syn ; Bartos: 0.09 ; Wang 0.1
para.ealpha_syn = 0.8; % Golomb & A
para.ebeta_syn = 0.5; %1/decay

%Others
para.eC_m = 1; %uF/cm^2
para.iC_m = 1; %uF/cm^2

%% colors for plotting
sett.red = [183 38 33]/256;
sett.blue = [56 100 176]/256;
sett.purple = [0.4 0 0.5];
sett.green = [0.1 0.6 0.2];

%% save settings;
sett.path = 'YourPathHere';
sett.saveloc = [sett.path, 'Figures/run',num2str(runnumber),'/run', num2str(runnumber), '_'];
sett.savelocdata = [sett.path, 'Data/run',num2str(runnumber),'/run', num2str(runnumber), '_'];

if savesettings == true
    mkdir([sett.path, 'Data/run',num2str(runnumber)]);
    mkdir([sett.path, 'Figures/run',num2str(runnumber)]);
    save([sett.savelocdata,'settings.mat'], 'sett','para');
    
    diary([sett.savelocdata, 'set.txt']);
    sett
    para
    diary off
    
    diary([sett.saveloc, 'set.txt']);
    sett
    para
    diary off
end

end