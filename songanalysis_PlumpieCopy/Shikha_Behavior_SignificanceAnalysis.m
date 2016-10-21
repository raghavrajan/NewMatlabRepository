function [] = Shikha_Behavior_SignificanceAnalysis(VocalizationTimesFile, FileDur, PreDur, PostDur, StimulusTimesFile)

% ========================================================================%
% First get the times of all vocalizations from the VocalizationTimesFile
% (text file)
% ========================================================================%

Fid = fopen(VocalizationTimesFile, 'r');
Temp = textscan(Fid, '%f', 'DeLimiter', '\n');
fclose(Fid);

VocalizationTimes = Temp{1};

% ========================================================================%
% First get the times of all stimuli from the StimulusTimesFile
% (text file)
% ========================================================================%

Fid = fopen(StimulusTimesFile, 'r');
Temp = textscan(Fid, '%s', 'DeLimiter', '\n');
fclose(Fid);

StimulusLines = find(cellfun(@length, strfind(Temp{1}, 'Stimulus')));
StimulusLines = [StimulusLines(:); (size(Temp{1},1) + 1)];
for i = 1:(length(StimulusLines) - 1),
    StimulusTimes{i} = str2double(Temp{1}((StimulusLines(i) + 1):(StimulusLines(i+1) - 1)));
end

NumStimuli = length(StimulusTimes{1});

% ========================================================================%
% Now choose random stimulus times and calculate the difference between
% number of vocalizations post-stimulus and pre-stimulus during a
% particular pre-stim period and post-stim period as specified by the
% inputs. Take the mean of this difference and repeat this procedure 1000
% times to get a distribution of mean differences that could arise purely
% by chance
% ========================================================================%

for i = 1:1000,
    RandomStimTimes = rand(NumStimuli, 1) * FileDur;
    clear Temp_PrePostDiff;
    for j = 1:length(RandomStimTimes),
        NumVoc_Pre = length(find((VocalizationTimes >= (RandomStimTimes(j) - PreDur)) & (VocalizationTimes < RandomStimTimes(j))));
        NumVoc_Post = length(find((VocalizationTimes >= RandomStimTimes(j)) & (VocalizationTimes < (RandomStimTimes(j) + PostDur))));
        Temp_PrePostDiff(j) = NumVoc_Post - NumVoc_Pre;
    end
    PrePostDiff(i) = mean(Temp_PrePostDiff);
end

PrePostDiff = sort(PrePostDiff);

for i = 1:length(StimulusTimes),
    clear Temp_PrePostDiff;
    for j = 1:length(StimulusTimes{i})
        NumVoc_Pre = length(find((VocalizationTimes >= (StimulusTimes{i}(j) - PreDur)) & (VocalizationTimes < StimulusTimes{i}(j))));
        NumVoc_Post = length(find((VocalizationTimes >= StimulusTimes{i}(j)) & (VocalizationTimes < (StimulusTimes{i}(j) + PostDur))));
        Temp_PrePostDiff(j) = NumVoc_Post - NumVoc_Pre;
    end
    ActualPrePostDiff{i} = mean(Temp_PrePostDiff);
end

figure;
set(gcf, 'Color', 'w', 'Position', [100 300 830 400]);
plot(linspace(min(PrePostDiff),max(PrePostDiff), 25), histc(PrePostDiff, linspace(min(PrePostDiff), max(PrePostDiff), 25)), 'k');
hold on;
axis tight;
temp = axis;

for i = 1:length(ActualPrePostDiff),
    plot([ActualPrePostDiff{i} ActualPrePostDiff{i}], temp(3:4), 'b--', 'LineWidth', 2);
    text(ActualPrePostDiff{i}, mean(temp(3:4)) * (1 + i*0.25), ['Actual difference for stim #', num2str(i), ' = ', num2str(ActualPrePostDiff{i})], 'Color', 'b', 'FontSize', 14);
end

Threshold = PrePostDiff(round(0.95*length(PrePostDiff)));
plot([Threshold Threshold], temp(3:4), 'r--');
text(Threshold, mean(temp(3:4)), ['95th percentile = ', num2str(Threshold)], 'Color', 'r', 'FontSize', 14);

LowerThreshold = PrePostDiff(round(0.05*length(PrePostDiff)));
plot([LowerThreshold LowerThreshold], temp(3:4), 'r--');
text(LowerThreshold, mean(temp(3:4)), ['5th percentile = ', num2str(LowerThreshold)], 'Color', 'r', 'FontSize', 14);

axis tight;
temp = axis;
temp(2) = 1.1*temp(2);
axis(temp);

xlabel('Difference between post-stimulus number of vocalizations and pre-stimulus number of vocalizations', 'FontSize', 14);
ylabel('#', 'FontSize', 14);
set(gca, 'FontSize', 14);
title('Distribution of differences that can arise with random stimulus times', 'FontSize', 16);

disp(['95% percentile of difference distribution is ', num2str(Threshold)]);

for i = 1:length(ActualPrePostDiff),
    disp(['Actual difference for stimulus #', num2str(i), ' is ', num2str(ActualPrePostDiff{i})]);
end


