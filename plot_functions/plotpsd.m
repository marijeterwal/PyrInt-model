function plotpsd(sett, f, spectrumi, spectrume, fstart, fend, e,i, saveind)

reg = sett.Nregions;

figure('Name','Spectra', 'Position',[100 100 600 400*(reg/1.5)]);

for r = 1:reg
subaxis(reg,2,1,r, 'SpacingVert',0.15,'SpacingHoriz', 0.10,'MarginTop',0.1,'MarginBottom',0.15, 'MarginRight',0.10,'MarginLeft',0.15, 'Padding',0);
plot(f,spectrume(:,r), 'Color', sett.red, 'LineWidth', 2);
xlim([fstart fend]);
xlabel ('Frequency (Hz)');
ylabel ({['Region ', num2str(r)], 'PSD (dB/Hz)'});
title('PSD of E');

subaxis(reg,2,2,r);
plot(f,spectrumi(:,r), 'Color', sett.blue, 'LineWidth', 2);
xlim([fstart fend]);
xlabel ('Frequency (Hz)');
ylabel ('PSD (dB/Hz)');
title('PSD of I');
end

if saveind == true
    saveas(gcf, [sett.saveloc,'spectra','_e', num2str(e),'_i', num2str(i)],'fig');
end