
% settings
if savefigures == true || savedata == true; savesettings = 1; else savesettings = 0; end

if ~exist('loaddata', 'var'); loaddata = 0; end
if ~exist('loadanalysis', 'var'); loadanalysis = 0; end

if loaddata == true || loadanalysis == true
    load([path, 'Data/run', num2str(runnumber),'/run', num2str(runnumber), '_','settings.mat']);
    sett.saveloc = [path, 'Figures/run',num2str(runnumber),'/run', num2str(runnumber), '_'];
    sett.savelocdata = [path, 'Data/run',num2str(runnumber),'/run', num2str(runnumber), '_'];
else
    [sett, para] = settings(savesettings,runnumber);
end

%% test settings

if sett.Tselection >= sett.Ttot; error('Tselection >= Ttot'); end

%% load common parameters from sett and para structures

% network parameters
dt = sett.dt;
Ttot = sett.Ttot;
tnumb = Ttot /dt +1;
tspan = 0:dt:Ttot;
Ni = sett.Ni;
Ne = sett.Ne;
reg = sett.Nregions;

% analysis parameters
Tselection = sett.Tselection;
corselection = sett.corselection;
dthist = sett.dthist;
cutoff_freq_low = sett.cutoff_freq_low;
cutoff_freq_high = sett.cutoff_freq_high;

% pulse pars
% ipulse = sett.Iipulse;
% epulse = sett.Iepulse;
% pulse_start = sett.Ipulse_start;
% pulse_amp = sett.Ipulse_amp;

% noise pars
if ~isfield(sett, 'noiseamp'); sett.noiseamp = 0.5; end
if ~isfield(sett, 'noisecor'); sett.noisecor = 0; end
% noiseamp_start = sett.noiseamp;

% periodic drive pars
if ~isfield(sett, 'Iper_freq'); sett.Iper_freq = 10; end
% Iper_period_start = sett.Iper_freq; %sett.Iper_period;
% Iper_amp_start = sett.Iper_amp;

if ~isfield(sett, 'g_syn_ifactor')
    sett.g_syn_ifactor = 1;
end

if ~isfield(sett, 'Inoise_loop')
    sett.Inoise_loop = 0;
    sett.noiseseed = sett.seed+ 19;
    sett.tau = 200; %in ms, for brownian noise
end
if ~isfield(sett, 'Idyn_frac'); sett.Idyn_frac = 0; end

% connection pars % only works for 2 areas
if length(sett.g_syn_ee_r(:,1)) == 1
    gsyn_ei_start = sett.g_syn_ei_r(1);
    gsyn_ee_start = sett.g_syn_ee_r(1);
    sett.g_syn_ee_r = zeros(reg,reg); sett.g_syn_ee_r(2,1) = gsyn_ee_start;
    sett.g_syn_ei_r = zeros(reg,reg); sett.g_syn_ei_r(2,1) = gsyn_ei_start;
else
    gsyn_ei_start = sett.g_syn_ei_r(2,1);
    gsyn_ee_start = sett.g_syn_ee_r(2,1);
end

if ~isfield(sett, 'Iestep')
    sett.Iistep = [0,0,0]; %
    sett.Iestep = [0,0,0]; %
end

if ~isfield(sett, 'g_syn_ii_r')
    sett.g_syn_ii_r = [0,0,0]; %
    sett.g_syn_ie_r = [0,0,0]; %
end


%% added parameters - move these to settings!

npeaks = 5;

params.Fs = 1000/dthist;
params.err = [2 0.05];
params.fpass = [sett.cutoff_freq_low, sett.cutoff_freq_high];
N = (Ttot-Tselection)/1000;
W = 30;
K = floor(2*N*W) - 1;
params.tapers = [N*W K];

nbinsMI = 6;

lateonset_starttime = sett.lateonset_start;

%% combinations of areas

if reg > 1
    comb = flipud(combnk(1:reg,2));
else
    comb = [1,1];
end
 ncomb = length(comb(:,1));