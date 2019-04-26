
%% synaptic delay
delayd = [3,4,5,6];
load('D:\Projects_CTC\fig5_syndel.mat') %syndel: seed x var x syn

%% frequency
freqd = [54,58,62,66,70,73];
load('D:\Projects_CTC\fig5_freqdel.mat') %syndel: seed x var x freq

%% E/I factor
fracd = [0.5, 0.75, 1, 1.25, 1.5];
load('D:\Projects_CTC\fig5_fracdel.mat') %syndel: seed x var x del

%% plot

orange = [232,78,27]/255;
yellow = [243,146,0]/255;

c1 = sett.purple; % firing rate
c2 = sett.green; % power

figure('Position', [50,50,1000,700])

% subplot 1
subplot(2,3,1); hold on % synaptic delay - phase
for i = 1:length(delayd)
    plot(delayd(i), syndel(:,2,i), 'o', 'color', c2, 'MarkerSize', 4);%, 'MarkerFaceColor', c2)
    plot(delayd(i), syndel(:,4,i), 'o', 'color', c1, 'MarkerSize', 4);%, 'MarkerFaceColor', c1);
end
p1 = polyfit(repmat(reshape(delayd,[1,1,4]),[5,1,1]), syndel(:,2,:), 1);
l2 = plot(delayd, p1(1)*delayd+p1(2), '-', 'color', c2, 'linewidth', 2);
p2 = polyfit(repmat(reshape(delayd,[1,1,4]),[5,1,1]), syndel(:,4,:), 1);
l1 = plot(delayd, p2(1)*delayd+p2(2), '-', 'color', c1, 'linewidth', 2);

datx = repmat(reshape(delayd,[1,1,4]),[5,1,1]);
daty = syndel(:,2,:);
[rho,p] = corr(datx(:), daty(:));
if p < 0.05
    plot(6,0.7, '*', 'color', c2, 'linewidth', 2);
end
datx = repmat(reshape(delayd,[1,1,4]),[5,1,1]);
daty = syndel(:,4,:);
[rho,p] = corr(datx(:), daty(:));
if p < 0.05
    plot(5.6,0.7, '*', 'color', c1, 'linewidth', 2);
end

xlim([2.5,6.5])
ylim([0.25,0.75])
ylabel('Optimal phase difference (rad)')
set(gca, 'ytick', [0.3,0.4,0.50,0.6,0.7] , 'yticklabel',{'0.6\pi', '0.8\pi','1\pi','1.2\pi','1.4\pi'})
set(gca, 'xticklabel', [])

% subplot 4
subplot(2,3,4); hold on % synaptic delay - amp
for i = 1:length(delayd)
    plot(delayd(i), syndel(:,1,i), 'o', 'color', c2, 'MarkerSize', 4);%, 'MarkerFaceColor', c2)
    plot(delayd(i), syndel(:,3,i), 'o', 'color', c1, 'MarkerSize', 4);%, 'MarkerFaceColor', c1);
end
p1 = polyfit(repmat(reshape(delayd,[1,1,4]),[5,1,1]), syndel(:,1,:), 1);
l2 = plot(delayd, p1(1)*delayd+p1(2), '-', 'color', c2, 'linewidth', 2);
p2 = polyfit(repmat(reshape(delayd,[1,1,4]),[5,1,1]), syndel(:,3,:), 1);
l1 = plot(delayd, p2(1)*delayd+p2(2), '-', 'color', c1, 'linewidth', 2);

datx = repmat(reshape(delayd,[1,1,4]),[5,1,1]);
daty = syndel(:,1,:);
[rho,p] = corr(datx(:), daty(:));
if p < 0.05
    plot(6,0.83, '*', 'color', c2, 'linewidth', 2);
end
datx = repmat(reshape(delayd,[1,1,4]),[5,1,1]);
daty = syndel(:,3,:);
[rho,p] = corr(datx(:), daty(:));
if p < 0.05
    plot(5.6,0.83, '*', 'color', c1, 'linewidth', 2);
end

legend([l1,l2],{'Firing rate', 'Power'}, 'Location', 'NorthWest')
legend boxoff
xlim([2.5,6.5])
ylim([0.07,0.9])
xlabel('Intercircuit synaptic delay (ms)')
ylabel('Max. correlation coefficient')

% subplot 2
subplot(2,3,2); hold on % frequency - phase
for i = 1:length(freqd)
    plot(freqd(i), freqdel(:,2,i), 'o', 'color', c2, 'MarkerSize', 4);%, 'MarkerFaceColor', c2)
    plot(freqd(i), freqdel(:,4,i), 'o', 'color', c1, 'MarkerSize', 4);%, 'MarkerFaceColor', c1);
end
p1 = polyfit(repmat(reshape(freqd,[1,1,6]),[5,1,1]), freqdel(:,2,:), 1);
plot(freqd, p1(1)*freqd+p1(2), '-', 'color', c2, 'linewidth', 2);
p2 = polyfit(repmat(reshape(freqd,[1,1,6]),[5,1,1]), freqdel(:,4,:), 1);
plot(freqd, p2(1)*freqd+p2(2), '-', 'color', c1, 'linewidth', 2);

datx = repmat(reshape(freqd,[1,1,6]),[5,1,1]);
daty = freqdel(:,2,:);
[rho,p] = corr(datx(:), daty(:))
if p < 0.05
    plot(72.5,0.7, '*', 'color', c2, 'linewidth', 2);
end
datx = repmat(reshape(freqd,[1,1,6]),[5,1,1]);
daty = freqdel(:,4,:);
[rho,p] = corr(datx(:), daty(:))
if p < 0.05
    plot(70,0.7, '*', 'color', c1, 'linewidth', 2);
end

xlim([50,75])
ylim([0.25,0.75])
set(gca, 'ytick', [0.3,0.4,0.50,0.6,0.7], 'yticklabel',{})
set(gca, 'xticklabel', [])

% subplot 5
subplot(2,3,5); hold on % frequency - amp
for i = 1:length(freqd)
    plot(freqd(i), freqdel(:,1,i), 'o', 'color', c2, 'MarkerSize', 4);%, 'MarkerFaceColor', c2)
    plot(freqd(i), freqdel(:,3,i), 'o', 'color', c1, 'MarkerSize', 4);%, 'MarkerFaceColor', c1);
end
p1 = polyfit(repmat(reshape(freqd,[1,1,6]),[5,1,1]), freqdel(:,1,:), 1);
l1 = plot(freqd, p1(1)*freqd+p1(2), '-', 'color', c2, 'linewidth', 2);
p2 = polyfit(repmat(reshape(freqd,[1,1,6]),[5,1,1]), freqdel(:,3,:), 1);
l2 = plot(freqd, p2(1)*freqd+p2(2), '-', 'color', c1, 'linewidth', 2);

datx = repmat(reshape(freqd,[1,1,6]),[5,1,1]);
daty = freqdel(:,1,:);
[rho,p] = corr(datx(:), daty(:));
if p < 0.05
    plot(72.5,0.83, '*', 'color', c2, 'linewidth', 2);
end
datx = repmat(reshape(freqd,[1,1,6]),[5,1,1]);
daty = freqdel(:,3,:);
[rho,p] = corr(datx(:), daty(:));
if p < 0.05
    plot(70,0.83, '*', 'color', c1, 'linewidth', 2);
end

xlim([50,75])
ylim([0.07,0.9])
xlabel('Oscillation frequency (Hz)')
set(gca, 'yticklabel', [])

% subplot 3
subplot(2,3,3); hold on % ei - phase
for i = 1:length(fracd)
    plot(fracd(i), fracdel(:,2,i), 'o', 'color', c2, 'MarkerSize', 4);%, 'MarkerFaceColor', c2)
    plot(fracd(i), fracdel(:,4,i), 'o', 'color', c1, 'MarkerSize', 4);%, 'MarkerFaceColor', c1);
end
p1 = polyfit(repmat(reshape(fracd,[1,1,5]),[5,1,1]), fracdel(:,2,:), 1);
plot(fracd, p1(1)*fracd+p1(2), '-', 'color', c2, 'linewidth', 2);
p2 = polyfit(repmat(reshape(fracd,[1,1,5]),[5,1,1]), fracdel(:,4,:), 1);
plot(fracd, p2(1)*fracd+p2(2), '-', 'color', c1, 'linewidth', 2);

datx = repmat(reshape(fracd,[1,1,5]),[5,1,1]);
daty = fracdel(:,2,:);
[rho,p] = corr(datx(:), daty(:));
if p < 0.05
    plot(1.5,0.7, '*', 'color', c2, 'linewidth', 2);
end
datx = repmat(reshape(fracd,[1,1,5]),[5,1,1]);
daty = fracdel(:,4,:);
[rho,p] = corr(datx(:), daty(:));
if p < 0.05
    plot(1.4,0.7, '*', 'color', c1, 'linewidth', 2);
end

xlim([0.4,1.6])
ylim([0.25,0.75])
set(gca, 'xtick', fracd, 'xticklabel', [])
set(gca, 'ytick', [0.3,0.4,0.50,0.6,0.7], 'yticklabel',{})

% subplot 6
subplot(2,3,6); hold on % ei - amp
for i = 1:length(fracd)
    plot(fracd(i), fracdel(:,1,i), 'o', 'color', c2, 'MarkerSize', 4);%, 'MarkerFaceColor', c2)
    plot(fracd(i), fracdel(:,3,i), 'o', 'color', c1, 'MarkerSize', 4);%, 'MarkerFaceColor', c1);
end
p1 = polyfit(repmat(reshape(fracd,[1,1,5]),[5,1,1]), fracdel(:,1,:), 1);
plot(fracd, p1(1)*fracd+p1(2), '-', 'color', c2, 'linewidth', 2);
p2 = polyfit(repmat(reshape(fracd,[1,1,5]),[5,1,1]), fracdel(:,3,:), 1);
plot(fracd, p2(1)*fracd+p2(2), '-', 'color', c1, 'linewidth', 2);

datx = repmat(reshape(fracd,[1,1,5]),[5,1,1]);
daty = fracdel(:,1,:);
[rho,p] = corr(datx(:), daty(:));
if p < 0.05
    plot(1.5,0.83, '*', 'color', c2, 'linewidth', 2);
end
datx = repmat(reshape(fracd,[1,1,5]),[5,1,1]);
daty = fracdel(:,3,:);
[rho,p] = corr(datx(:), daty(:));
if p < 0.05
    plot(1.4,0.83, '*', 'color', c1, 'linewidth', 2);
end

xlim([0.4,1.6])
ylim([0.07,0.9])
set(gca, 'xtick', fracd)
set(gca, 'yticklabel', [])
xlabel('E\rightarrowI / E\rightarrowE cond. ratio')

