function plotautocorrelation(set, auto_i, auto_e, tspan,e, i, saveind)

reg = set.Nregions;

figure('Name','Correlations', 'Position',[100 100 600 400*(reg/1.5)]);

for r = 1:reg
    subaxis(reg,2,1,r, 'SpacingVert',0.15,'SpacingHoriz', 0.10,'MarginTop',0.1,'MarginBottom',0.15, 'MarginRight',0.10,'MarginLeft',0.15, 'Padding',0);
    h1 = area(tspan, auto_e(:,:,r), 'FaceColor', set.red, 'EdgeColor', set.red);
    ylim([-1 1]);
    title(['Region', num2str(r)]);
    if r == 1; ylabel ('Autocorrelation of E'); end
        
    subaxis(reg,2,2,r);
    h3 = area(tspan, auto_i(:,:,r), 'FaceColor', set.blue, 'EdgeColor', set.blue);
    if r == 1; ylabel ('Autocorrelation of I'); end;
    ylim([-1 1]);
    xlabel ('Delay (ms)');
end

if saveind == true
    saveas(gcf, [set.saveloc,'correlations','_e', num2str(e),'_i', num2str(i)],'fig');
end

end