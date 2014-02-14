function [idx] = PlotFeatures(SpikeTimes, Waveforms, FileName, ClusterNos, PeakIndex, Fs, UpSampleFactor, LogFid)

disp(['Peak index is ', num2str(PeakIndex)]);

Wv1 = reshape(Waveforms(:,1,:), size(Waveforms,1), size(Waveforms,3));
Wv1 = Wv1';

PeakInds = (PeakIndex - 5*UpSampleFactor):(PeakIndex + 3*UpSampleFactor);
PreVallInds = 1:(PeakIndex - 5*UpSampleFactor);
if (size(Waveforms,3) > 32)
    PostVallInds = (PeakIndex + 3*UpSampleFactor):(PeakIndex + 18*UpSampleFactor);
else
    PostVallInds = (PeakIndex + 3*UpSampleFactor):(size(Waveforms,3) - 2);
end
RisSlopeInds = 1:PeakIndex;
FallSlopeInds = PeakIndex:PostVallInds(end);

peak = max(Wv1(PeakInds,:))';
[tempPeak, peakIndex] = max(Wv1(PeakInds,:));
pre_valley = min(Wv1(PreVallInds,:))';
post_valley = min(Wv1(PostVallInds,:))';
[temppost_valley, postValleyIndex] = min(Wv1(PostVallInds,:));
energy = (sqrt(sum(Wv1.*Wv1))/size(Wv1,1))';
rising_slope = max(diff(Wv1(RisSlopeInds,:)))';
falling_slope = min(diff(Wv1(FallSlopeInds,:)))';
spike_width = postValleyIndex - peakIndex;
if (size(spike_width,1) < size(spike_width,2))
    spike_width = spike_width';
end


[coeff, score, latent] = princomp(Wv1');
for i = 1:10,
    disp(['PC', num2str(i), ' and its contribution is ', num2str((latent(i)/sum(latent))*100)]);
    fprintf(LogFid, 'PC%i and its contribution is %g\n', i, (latent(i)/sum(latent))*100);
end

FeatureVect = [score(:,1:5) peak];
fprintf(LogFid, 'Features used for clustering are PC1 to PC5 and peak\n');

[idx] = kmeans(FeatureVect, ClusterNos, 'replicates', 200);

Markers = ['b+'; 'r+'; 'g+'; 'k+'; 'c+'; 'm+'; 'y+'];
Colors{1} = 'blue';
Colors{2} = 'red';
Colors{3} = 'green';
Colors{4} = 'black';
Colors{5} = 'cyan';
Colors{6} = 'magenta';
Colors{7} = 'yellow';

figure(1);
subplot(3,3,1);
hold off;
for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(peak(idx==i), post_valley(idx==i), Markers(mod(i,7),:), 'MarkerSize', 2);
    hold on;
end
title(['Peak vs. Post Valley '], 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

subplot(3,3,2);
hold off;
for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(score(idx==i, 1), score(idx==i, 2), Markers(mod(i,7),:), 'MarkerSize', 2);
    hold on;
end
title(['PC1 vs. PC2'], 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

subplot(3,3,3);
hold off;
for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(peak(idx==i), score(idx==i, 1), Markers(mod(i,7),:), 'MarkerSize', 2);
    hold on;
end
title(['Pk vs PC1 '], 'FontSize', 12, 'FontWeight', 'bold');
axis tight;
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

subplot(3,3,4);
hold off;
for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(post_valley(idx==i), score(idx==i, 1), Markers(mod(i,7),:), 'MarkerSize', 2);
    hold on;
end
title(['Post Valley vs PC1 '], 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

subplot(3,3,5);
hold off;
for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(peak(idx==i), rising_slope(idx==i), Markers(mod(i,7),:), 'MarkerSize', 2);
    hold on;
end
title(['Pk vs Max rising slope'], 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

subplot(3,3,6);
hold off;
for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(peak(idx==i), energy(idx==i), Markers(mod(i,7),:), 'MarkerSize', 2);
    hold on;
end
title(['Pk vs Energy'], 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

subplot(3,3,7);
hold off;
for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(score(idx==i, 1), energy(idx==i), Markers(mod(i,7),:), 'MarkerSize', 2);
    hold on;
end
title(['PC1 vs Energy'], 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

subplot(3,3,8);
hold off;
for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(peak(idx==i), pre_valley(idx==i), Markers(mod(i,7),:), 'MarkerSize', 2);
    hold on;
end
title(['Pk vs Pre Valley'], 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

subplot(3,3,9);
hold off;
for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(peak(idx==i), falling_slope(idx==i), Markers(mod(i,7),:), 'MarkerSize', 2);
    hold on;
end
title(['Pk vs Falling Slope'], 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
set(gcf, 'Color', 'w');

figure(2);
WaveformAxis = axes('Position', [0.1 0.12 0.4 0.8]);
hold off;
WaveformIndex = 1000*(1:1:size(Wv1, 1))/Fs;

for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(WaveformIndex, Wv1(:,idx==i)*100, Markers(mod(i,7),1));
    hold on;
end
title(['Waveforms ', FileName], 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Peak spike amplitude (uV)', 'FontSize', 14, 'FontWeight', 'bold');
set(gcf, 'Color', 'w');

PlotHeight = 0.85/(max(idx) * 1.1);

for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    
    ISIAxis = axes('Position', [0.55 (0.12 + ((i-1)*PlotHeight) +((i-1)*0.05)) 0.4 PlotHeight]);
    Edges = 0:0.1:1000;
    ISI = histc(diff(SpikeTimes(idx==i) * 1000), Edges);
    ISI = ISI/sum(ISI);
    ISIPlot = bar(Edges, ISI*100);
    set(ISIPlot, 'FaceColor', 'w', 'EdgeColor', Markers(mod(i,7),1));
    set(ISIAxis, 'XScale', 'log');
    set(ISIAxis, 'FontSize', 12, 'FontWeight', 'bold');
    axis([0 1000 0 (max(ISI)*110)]);
end

figure(3);
hold off;
for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(SpikeTimes(idx == i), peak(idx == i)*100, Markers(mod(i,7),:));
    hold on;
end
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time (sec)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Peak spike amplitude (uV)', 'FontSize', 12, 'FontWeight', 'bold');
set(gcf, 'Color', 'w');
title('Peak amplitude vs time', 'FontSize', 12, 'FontWeight', 'bold');

for i = 1:max(idx),
    Results{i} = ['No of points in ', Colors{i}, ' cluster is ', num2str(length(find(idx == i)))];
end

for i = 1:max(idx),
    try
        IsolDist(i) = IsolationDistance(FeatureVect, find(idx == i));
    catch
        IsolDist(i) = -100;
    end
    try
        [L, Lr, df] = L_Ratio(FeatureVect, find(idx == i));
        Lratio(i) = Lr;
    catch
        Lratio(i) = -100;
    end
end

for i = 1:max(idx),
    Results{max(idx) + i} = [Colors{i}, ' cluster: Isolation distance is ', num2str(IsolDist(i)), ' and L-ratio is ', num2str(Lratio(i))];
end

msgbox(Results, 'Results');
