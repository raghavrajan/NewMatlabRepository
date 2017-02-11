function [idx] = PlotTetrodeFeatures(SpikeTimes, Waveforms, FileName, ClusterNos, PeakIndex, Fs)

Wv1 = reshape(Waveforms(:,1,:), size(Waveforms,1), size(Waveforms,3));
Wv1 = Wv1';
Wv2 = reshape(Waveforms(:,2,:), size(Waveforms,1), size(Waveforms,3));
Wv2 = Wv2';
Wv3 = reshape(Waveforms(:,3,:), size(Waveforms,1), size(Waveforms,3));
Wv3 = Wv3';
Wv4 = reshape(Waveforms(:,4,:), size(Waveforms,1), size(Waveforms,3));
Wv4 = Wv4';

PeakInds = (PeakIndex - 5):(PeakIndex + 3);
if (size(Waveforms,1) > 32)
    PostVallInds = (PeakIndex + 3):(PeakIndex + 18);
else
    PostVallInds = (PeakIndex + 3):(size(Waveforms,1) - 2);
end
RisSlopeInds = 1:PeakIndex;
FallSlopeInds = PeakIndex:PostVallInds(end);

peak1 = max(Wv1(PeakInds,:))';
peak2 = max(Wv2(PeakInds,:))';
peak3 = max(Wv3(PeakInds,:))';
peak4 = max(Wv4(PeakInds,:))';

%energy = (sqrt(sum(Wv1.*Wv1))/size(Wv1,1))';
%rising_slope = max(diff(Wv1(RisSlopeInds,:)))';
%falling_slope = min(diff(Wv1(FallSlopeInds,:)))';

[coeff, score, latent] = princomp(Wv1');

PCScores(:,1) = Wv1'*coeff(:,1);
PCScores(:,2) = Wv2'*coeff(:,1);
PCScores(:,3) = Wv3'*coeff(:,1);
PCScores(:,4) = Wv4'*coeff(:,1);

FeatureVect = [PCScores peak1 peak2 peak3 peak4];

[idx] = kmeans(FeatureVect, ClusterNos, 'replicates', 1000);

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
    plot(FeatureVect(idx==i,1), FeatureVect(idx==i,2), Markers(mod(i,7),:), 'MarkerSize', 2);
    hold on;
end
title(['PC1 vs. PC2 '], 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

subplot(3,3,2);
hold off;
for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(FeatureVect(idx==i, 1), FeatureVect(idx==i, 3), Markers(mod(i,7),:), 'MarkerSize', 2);
    hold on;
end
title(['PC1 vs. PC3'], 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
 
subplot(3,3,3);
hold off;
for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(FeatureVect(idx==i,1), FeatureVect(idx==i, 4), Markers(mod(i,7),:), 'MarkerSize', 2);
    hold on;
end
title(['PC1 vs PC4 '], 'FontSize', 12, 'FontWeight', 'bold');
axis tight;
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

subplot(3,3,4);
hold off;
for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(FeatureVect(idx==i,2), FeatureVect(idx==i, 3), Markers(mod(i,7),:), 'MarkerSize', 2);
    hold on;
end
title(['PC2 vs PC3 '], 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

subplot(3,3,5);
hold off;
for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(FeatureVect(idx==i,2), FeatureVect(idx==i,4), Markers(mod(i,7),:), 'MarkerSize', 2);
    hold on;
end
title(['PC2 vs PC4'], 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

subplot(3,3,6);
hold off;
for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(FeatureVect(idx==i,3), FeatureVect(idx==i,4), Markers(mod(i,7),:), 'MarkerSize', 2);
    hold on;
end
title(['PC3 vs PC4'], 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

subplot(3,3,7);
hold off;
for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(FeatureVect(idx==i, 1), FeatureVect(idx==i,5), Markers(mod(i,7),:), 'MarkerSize', 2);
    hold on;
end
title(['PC1 vs Peak1'], 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

subplot(3,3,8);
hold off;
for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(FeatureVect(idx==i,1), FeatureVect(idx==i,3), Markers(mod(i,7),:), 'MarkerSize', 2);
    hold on;
end
title(['PC1 vs Peak2'], 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

subplot(3,3,9);
hold off;
for i = 1:max(idx),
    if (isempty(find(idx==i, 1)))
        continue;
    end
    plot(FeatureVect(idx==i,1), FeatureVect(idx==i,3), Markers(mod(i,7),:), 'MarkerSize', 2);
    hold on;
end
title(['PC1 vs Peak3'], 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
set(gcf, 'Color', 'w');

% figure(2);
% WaveformAxis = axes('Position', [0.1 0.12 0.4 0.8]);
% hold off;
% WaveformIndex = 1000*(1:1:size(Wv1, 1))/Fs;
% 
% for i = 1:max(idx),
%     if (isempty(find(idx==i, 1)))
%         continue;
%     end
%     plot(WaveformIndex, Wv1(:,idx==i)*100, Markers(mod(i,7),1));
%     hold on;
% end
% title(['Waveforms ', FileName], 'FontSize', 12, 'FontWeight', 'bold');
% set(gca, 'FontSize', 12, 'FontWeight', 'bold');
% xlabel('Time (ms)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Peak spike amplitude (uV)', 'FontSize', 12, 'FontWeight', 'bold');
% set(gcf, 'Color', 'w');
% 
% PlotHeight = 0.85/(max(idx) * 1.1);
% 
% for i = 1:max(idx),
%     if (isempty(find(idx==i, 1)))
%         continue;
%     end
%     
%     ISIAxis = axes('Position', [0.55 (0.12 + ((i-1)*PlotHeight) +((i-1)*0.05)) 0.4 PlotHeight]);
%     Edges = 0:0.1:1000;
%     ISI = histc(diff(SpikeTimes(idx==i) * 1000), Edges);
%     ISI = ISI/sum(ISI);
%     bar(Edges, ISI);
%     set(gca, 'XScale', 'log');
%     set(gca, 'FontSize', 12, 'FontWeight', 'bold');
%     axis([0 1000 0 (max(ISI)*1.1)]);
% end
% 
% figure(3);
% hold off;
% for i = 1:max(idx),
%     if (isempty(find(idx==i, 1)))
%         continue;
%     end
%     plot(SpikeTimes(idx == i), peak(idx == i)*100, Markers(mod(i,7),:));
%     hold on;
% end
% set(gca, 'FontSize', 12, 'FontWeight', 'bold');
% xlabel('Time (sec)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Peak spike amplitude (uV)', 'FontSize', 12, 'FontWeight', 'bold');
% set(gcf, 'Color', 'w');
% title('Peak amplitude vs time', 'FontSize', 12, 'FontWeight', 'bold');
% 
% for i = 1:max(idx),
%     Results{i} = ['No of points in ', Colors{i}, ' cluster is ', num2str(length(find(idx == i)))];
% end
% 
% for i = 1:max(idx),
%     try
%         IsolDist(i) = IsolationDistance(FeatureVect, find(idx == i));
%     catch
%         IsolDist(i) = -100;
%     end
%     try
%         [L, Lr, df] = L_Ratio(FeatureVect, find(idx == i));
%         Lratio(i) = Lr;
%     catch
%         Lratio(i) = -100;
%     end
% end
% 
% for i = 1:max(idx),
%     Results{max(idx) + i} = [Colors{i}, ' cluster: Isolation distance is ', num2str(IsolDist(i)), ' and L-ratio is ', num2str(Lratio(i))];
% end
% 
% msgbox(Results, 'Results');
