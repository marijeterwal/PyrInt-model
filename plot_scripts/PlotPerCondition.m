% plots per condition

%% settings
plotraster = 1;
plotspikehistogram = 1;
plotautocor = 0;
plotcrosscor = 0;
plotspectra = 0;
plotcoherenceandphase = 0;
plotpulse = 0;

plotselection_start = 400; %sett.Tselection + (sett.windowsize/2);%350; % ms
plotselection_end = 700; %sett.Ttot - (sett.windowsize/2);%600; % ms

%% plots

if plotraster == true; plotrastergram(sett, spikes_i, spikes_e, plotselection_start, plotselection_end, l1,l2, savefigures); end

if plotspikehistogram == true; t_start = plotselection_start/dthist; t_end = plotselection_end/dthist;
     plotsth(sett, hist_i, hist_e, histtime_i,plotselection_start, plotselection_end,  l1,l2, savefigures);
end

if plotautocor == true;
    if timeresolved == false;
        plotautocorrelation(sett, corii, coree, 0:dthist:corselection,  l1,l2, savefigures);
    else plotautocorrelationtimeresolved(sett, corii_timeres, 0:dthist:corselection, windowtimes, slopeii_timeres, avgcorii_timeres,  l1,l2, savefigures);
    end
end

if plotcrosscor == true;
    cortimes = -(sett.Ttot - Tselection):dthist:(sett.Ttot - Tselection);
    if timeresolved == false;
        if reg == 1; plotcrosscorrelation(sett, corei, cortimes,  l1,l2, 'Correlation between E and I', reg, savefigures);
        elseif reg == 2; plotcrosscorrelation(sett, cor12, cortimes,  l1,l2, 'Correlation between region 1 and 2',  1, savefigures);
        end
    elseif timeresolved == true && reg == 2;
        cortimes = -(windowsize * dthist) : dthist:(windowsize * dthist);
        plotcrosscorrelationtimeresolved(sett, cor12_timeres, cortimes, windowtimes, slope12_timeres, avgcor12_timeres,  l1,l2, savefigures);
    end
end

if plotspectra == true;
    if timeresolved == false; plotpsd(sett, fpsd, psd_i, psd_e, cutoff_freq_low, cutoff_freq_high,  l1,l2, savefigures);
    else plotpsdtimeresolved(sett, psdi_timeres, fpsd_timeres, windowtimes, cutoff_freq_low, cutoff_freq_high, avg_psdi_timeres, peakfreq_psdi_timeres, width_psdi_timeres, e, i, savefigures);
    end
end

if plotcoherenceandphase == true;
    if timeresolved == false; plotcoherence(sett, cohe, phase, fcoh, cutoff_freq_low, cutoff_freq_high,  l1,l2, savefigures);
        %                 plotcoherence(sett, cohold, phaseold, fcohold, cutoff_freq_low, cutoff_freq_high, e, i, savefigures);
    else plotcoherencetimeresolved(sett, coh_timeres, phase_timeres, fcoh_timeres, windowtimes, cutoff_freq_low, cutoff_freq_high,  l1,l2, savefigures);
    end
end

if plotpulse == true
    plotstart = (sett.Ipulse_start - 30)/dthist;
    plotend = (sett.Ipulse_start + 70)/dthist;

    figure(100)
    subplot(nsubplots, 4, (nplot-1)*4+1); hold on
    plot(histtime_e(1,plotstart:plotend,1), hist_e0(1,plotstart:plotend,1), 'color', [0.6,0.6,0.6]);
    plot(histtime_e(1,plotstart:plotend,1), hist_e(1,plotstart:plotend,1), 'color', sett.red);
    plot([pulse_start, pulse_start], [0, 400], 'k')
    xlim([plotstart, plotend]*dthist);
    ylabel(['\theta = ', num2str(phaseOfPulse(l1,l2,1), 2)]);
    
    subplot(nsubplots, 4, (nplot-1)*4+2); hold on
    plot(histtime_e(1,plotstart:plotend,1), freqosc0(1,plotstart-20:plotend-20,1), 'color', [0.6,0.6,0.6]);
    plot(histtime_e(1,plotstart:plotend,1), freqosc(1,plotstart-20:plotend-20,1), 'color', sett.blue);
    plot([pulse_start, pulse_start], [70, 80], 'k')
    xlim([plotstart, plotend]*dthist);

    subplot(nsubplots, 4, (nplot-1)*4+3); hold on
    plot(histtime_e(1,plotstart:plotend,1), sigpeak0(1,plotstart-20:plotend-20,1)-410, 'color', [0.6,0.6,0.6]);
    plot(histtime_e(1,plotstart:plotend,1), sigpeak(1,plotstart-20:plotend-20,1)-410, 'color', sett.green);
    plot([pulse_start, pulse_start], [000, 400], 'k')
    xlim([plotstart, plotend]*dthist);

    subplot(nsubplots, 4, (nplot-1)*4+4); hold on
    plot(histtime_e(1,plotstart:plotend,1), sigbase0(1,plotstart-20:plotend-20,1), 'color', [0.6,0.6,0.6]);
    plot(histtime_e(1,plotstart:plotend,1), sigbase(1,plotstart-20:plotend-20,1), 'color', sett.purple);
    plot([pulse_start, pulse_start], [70, 90], 'k')
    xlim([plotstart, plotend]*dthist);
end
