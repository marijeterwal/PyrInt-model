function plotsth_regions(sett, hist_i, hist_e, tspan, tstart, tend, e,i, saveind)

reg = sett.Nregions;

figure('Name','STH', 'Position',[100 100 600 400*(reg/1.5)]);

for r = 1:reg
    subaxis(reg,1,r, 'SpacingVert',0.15,'SpacingHoriz', 0.05,'MarginTop',0.1,'MarginBottom',0.15, 'MarginRight',0.10,'MarginLeft',0.2, 'Padding',0);
    h1 = area(tspan, hist_e(:,:,r), 'FaceColor', sett.red, 'EdgeColor', sett.red);
    hold on;
    h3 = area(tspan, -hist_i(:,:,r), 'FaceColor', sett.blue, 'EdgeColor', sett.blue);
    if r == 1; legend ([h1,h3], {'Pyr.','Int.'},'Location','NorthEast'); end
    hold off;
    ylabel ({['Spike rate (spikes/s)']});
    ylim ([-250 250]);
    set(gca, 'Ytick', [-200 -100 0 100 200], 'YTickLabel','200|100|0|100|200');
    xlim ([tstart tend]);
%     title(['Region ', num2str(r)]);
end

xlabel ('Time (ms)');

if saveind == true
    saveas(gcf, [sett.saveloc,'STH','_e', num2str(e),'_i', num2str(i)],'fig');
    %     saveas(gcf, [set.saveloc,num2str(i),'_smoothed_spike_time_hist'],'pdf');
end

end