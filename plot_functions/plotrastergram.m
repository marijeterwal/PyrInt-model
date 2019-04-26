function plotrastergram(sett, Yi, Ye, tstart, tend, e, i, saveind)
% requires the datasets to be: [neuron, times,reg]
N = sett.Ni + sett.Ne;
reg = sett.Nregions;
figure('Name','Rastergram', 'Position',[100 100 600 500*(reg/1.5)]);
hold on;

% if reg == 1
% scatter(Yi(:,2),Yi(:,1),'Marker', '.', 'MarkerEdgeColor', set.blue, 'MarkerFaceColor', set.blue);
% hold on;
% scatter(Ye(:,2),set.Ni + Ye(:,1),'Marker', '.', 'MarkerEdgeColor', set.red, 'MarkerFaceColor', set.red); %
% xlim ([tstart tend]);
% ylim ([0 N]);
% ylabel ('Neurons');
% xlabel ('Time (ms)');
% legend ('Inh', 'Exc','Location','NorthEast');
% else
    for r = 1:reg
        subplot(reg,1,r); hold on;
        for n = 1:sett.Ni
            plot([Yi(Yi(:,3) == r & Yi(:,1) == n,2), Yi(Yi(:,3) == r & Yi(:,1) == n,2)]', [Yi(Yi(:,3) == r& Yi(:,1) == n,1)-0.5, Yi(Yi(:,3) == r & Yi(:,1) == n,1)+0.5]', 'color', sett.blue, 'linewidth', 1.5);
        end
        for n = 1:sett.Ne
            plot([Ye(Ye(:,3) == r & Ye(:,1) == n,2), Ye(Ye(:,3) == r & Ye(:,1) == n,2)]', [sett.Ni + Ye(Ye(:,3) == r & Ye(:,1) == n,1)-0.5, sett.Ni + Ye(Ye(:,3) == r & Ye(:,1) == n,1)+0.5]', 'color', sett.red, 'linewidth', 1.5);
        end
%         scatter(Yi(Yi(:,3) == r,2),Yi(Yi(:,3) == r,1),'Marker', '.', 'MarkerEdgeColor', sett.blue, 'MarkerFaceColor', sett.blue);
%         scatter(Ye(Ye(:,3) == r,2),sett.Ni + Ye(Ye(:,3) == r,1),'Marker', '.', 'MarkerEdgeColor', sett.red, 'MarkerFaceColor', sett.red);
        ylabel ('Neuron #');
        xlim ([tstart tend]);
        ylim ([0 N]);
        title(['Region ', num2str(r)]);
    end
    xlabel ('time (ms)');
% end

if saveind == true
        saveas(gcf, [sett.saveloc,'rastergram','_e', num2str(e),'_i', num2str(i)],'fig');
end
end