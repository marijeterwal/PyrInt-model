function savespikes(saveloc, spikes_i, spikes_e, l1,l2)
save([saveloc, 'l1', num2str(l1),'_l2', num2str(l2), '.mat'], 'spikes_i', 'spikes_e');
end