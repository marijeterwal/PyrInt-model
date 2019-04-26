function plotcrosscorrelation(set, cor, tspan, e, i, titlestr, reg, saveind)

figure('Name','Crosscorrelations', 'Position',[100 100 600 400*(4/3)]);

for r = 1:reg
    subaxis(1,reg,r, 'SpacingVert',0.15,'SpacingHoriz', 0.10,'MarginTop',0.1,'MarginBottom',0.15, 'MarginRight',0.10,'MarginLeft',0.15, 'Padding',0);
    h1 = area(tspan, cor(:,:,r), 'FaceColor', set.purple, 'EdgeColor', set.purple);
    ylim([-1 1]);
    xlabel('Delay (ms)');
%     title(['Region', num2str(r)]);
title(titlestr);
    if r == 1; ylabel ('Crosscorrelation'); end
end

if saveind == true
    saveas(gcf, [set.saveloc,'crosscorrelations','_e', num2str(e),'_i', num2str(i)],'fig');
end

end