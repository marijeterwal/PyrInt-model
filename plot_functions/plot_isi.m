function plot_isi(sett, spikes_i, spikes_e, e,i, saveind)

reg = sett.Nregions;

spikes_itemp = spikes_i;
spikes_itemp(:,1) = spikes_itemp(:,1) + sett.Ne;
N = Ne+Ni;

spikes = [spikes_itemp; spikes_e];
% spikes = spikes_e;
binc = linspace(0,100,200); %ms
count = zeros(N,length(binc));

for n = 1:N;
    [count(n,:), xout] = hist(diff(spikes(spikes(:,1) == n, 2)), binc);
end

count = sum(count, 1);

figure('Name','ISI', 'Position',[100 100 600 400*(reg/1.5)]);
subaxis(1,1,1, 'SpacingVert',0.15,'SpacingHoriz', 0.05,'MarginTop',0.1,'MarginBottom',0.15, 'MarginRight',0.10,'MarginLeft',0.2, 'Padding',0);
bar(binc, count, 1, 'FaceColor', sett.purple, 'Edgecolor', sett.purple);
xlim([0 60]);
ylim([0 5000]);
ylabel('Count');
xlabel('Interspike Interval (ms)');

saveas(gcf, [sett.saveloc,'ISI'],'fig');
saveas(gcf, [sett.saveloc, 'ISI'],'pdf');

end

