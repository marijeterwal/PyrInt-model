function plotcoherence(set, coh, phase, F, fstart, fend, e, i, saveind)

if fstart < F(2); fs = F(2);
else fs = fstart;
end

figure('Name','Coherence', 'Position',[100 100 600 400*(4/3)]);

    subaxis(2,1,1, 'SpacingVert',0.15,'SpacingHoriz', 0.05,'MarginTop',0.1,'MarginBottom',0.15, 'MarginRight',0.10,'MarginLeft',0.15, 'Padding',0);
    plot(F(2:end), coh(2:end), 'Color', set.purple , 'LineWidth', 2);%, 'LineWidth', 2); %fspan(2:end),
        ylabel ('Coherence');
        xlim([fs fend]);
        
    subaxis(2,1,2);
    plot(F(2:end), phase(2:end), 'Color', set.green ,'LineWidth', 2); %fspan(2:end), 
    ylabel ('Phase (2\pi rad)');
     xlim([fs fend]);


xlabel ('Frequency (Hz)');

if saveind == true
    saveas(gcf, [set.saveloc,'coherence','_e', num2str(e),'_i', num2str(i)],'fig');
end

end