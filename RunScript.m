% RunScript

% The data generation can be run in parallel, using the parfor for the l1 loop ...
% or in the case of applied pulses, with parfor for the l2 loop.
% The analysis cannot be done in parallel.

%% preliminary...

clear all
clc
close all

%% Loading and saving - settings

loaddata = 0; % 0 = calc data; 1 = load from file
loadanalysis = 0; % 0 = redo analysis; 1 = load from file
% if load, set path:
path = 'YourPathHere';

runnumber = 2;

savefigures = 0;
savedata = 1;


%% Analyses and plots - settings

analyzePPC = 1;
analyzeTimeRes   = 0;
analyzeInfo      = 0;

plotDrive = 1; % measures vs input current for one area, vs freqdiff for two areas
plotPerCondition = 0; % see script for further specification

%% parameters and settings  - initialization

GetPars % load common parameters
% sett.Iapploop = 0;
sett_ini = sett;

DetermineLoops % determine which parameters are looped over
Preallocation % allocate memory for analysis

%% load some currents
LoadNoise

if sum(sett.Iipulse + sett.Iepulse) > 0; analyzePulses = 1;
else analyzePulses = 0;
end

%% Run network

if loadanalysis == true
    fprintf('Loading analyses...');
    LoadAnalysis
    fprintf(' Done\n');
else
    
    for l1 = 1:l1steps
        
        if analyzePulses == true && sum(sett.Ieper+sett.Iiper)==0
            l2 = 1;
            csett = changePars(sett, loopId, l1, l2); % change l1 parameter
            ReferenceRun;
        end
        
        for l2 = 1:l2steps
            
            if analyzePulses == true && sum(sett.Ieper+sett.Iiper)
                csett = changePars(sett, loopId, l1, l2); % change l2 parameter
                ReferenceRun;
            end
            
            csett = changePars(sett, loopId, l1, l2); % change l2 parameter
            
            % load or run
            tic;
            if loaddata == true % Load data
                fprintf('Loading data: l1 = %i of %i; l2 = %i of %i\n', l1, l1steps, l2, l2steps);
                data = open([sett.savelocdata, 'l1', num2str(l1),'_l2', num2str(l2),'.mat']); %run 660
                spikes_i = data.spikes_i; spikes_e = data.spikes_e;
            else % Run network
                fprintf('Calculating data: l1 = %i of %i; l2 = %i of %i\n', l1, l1steps, l2, l2steps);
                [spikes_i, spikes_e] = calcEInetworkRK4(sett, para, csett);
                if savedata == true; savespikes(sett.savelocdata, spikes_i, spikes_e, l1, l2); end
            end
            toc
            
            % analysis - not performed when using Parfor
            AnalysisStandard1
            if reg > 1;                  AnalysisStandard2;      end
            if analyzePulses == true && length(pulse_start) == 1
                AnalysisPulses;                                  end
            if analyzeInfo == true;      AnalysisInfo;           end
            
            % plots
            if plotPerCondition == true; PlotPerCondition;       end
        end
    end
    
    %% Save
    SaveAnalysis
    
end

%% Plot

if plotDrive == true
    if reg == 1
        plotFreqe = freqe';
        plotFreqi = freqi';
        meanRe = re';
        meanRi = ri';
        plotPle = Phaselocke';
        plotPli = Phaselocki';
        plotCV   = CVe';
        PlotDriveOneArea
    else
        variableFreqAxis = 0;
        meanPli_12 = Phaselocki_12;
        meanPle_12 = Phaselocke_12;
        Phaseli = Phaselocki;
        Phasele = Phaselocke;
        PlotDriveTwoAreas
    end
end

%Pulses
if analyzePulses == true
    PlotPulseOneArea
    if reg > 1
        PlotPulseTwoAreas
    end
end

% Information
if analyzeInfo == true && reg == 2 && l1steps > 1
    PlotInfoTwoAreas
end
