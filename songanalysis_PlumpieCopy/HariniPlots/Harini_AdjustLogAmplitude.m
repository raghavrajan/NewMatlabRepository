function [RawSyllLogAmplitudeMeanValue, AdjustedSyllLogAmplitudeMeanValue] = Harini_AdjustLogAmplitude(IndividualBirds, BirdNames)

% This function is to calculate aronov fee amplitudes for all syllables and
% then adjust them across days to a common threshold based on the max of
% the minimum pre-data across days/sessions

ThresholdPercentile = 99;
PrePostPadding = 0.01; % in seconds
FHigh = 860;
FLow = 8600;
FFTWinSize = 5; % in ms

for Bird = 1:length(IndividualBirds),
    TempRawSyllLogAmplitudes = [];
    TempAdjustedSyllLogAmplitudes = [];
    fprintf('\n%s :\n', BirdNames);
    IndividualBirds(Bird).SyllLogAmplitude = [];
    IndividualBirds(Bird).RawSyllLogAmplitudeMeanValue = [];
    IndividualBirds(Bird).AdjustedSyllLogAmplitudeMeanValue = [];
    IndividualBirds(Bird).SyllLogAmplitudeStatus = [];
    FilesWithSylls = unique(IndividualBirds(Bird).AllSyllableData(:,2));
    for i = FilesWithSylls(:)',
        fprintf('%d > ', i);
        [DataDir, FileName, Ext] = fileparts(IndividualBirds(Bird).SongFileNames{i});
        [RawData, Fs] = GetData(DataDir, [FileName, Ext], IndividualBirds(Bird).SortedBirdParameters(1).FileType, 0);
        
        LogAmplitude = ASSLCalculateLogAmplitudeAronovFee(RawData, Fs, [], FFTWinSize, [], FHigh, FLow);

        % Now all sylls in that file
        Sylls = find(IndividualBirds(Bird).AllSyllableData(:,2) == i);
        % First Raw Waveform
        for j = Sylls(:)',
            StartIndex = round(((IndividualBirds(Bird).AllSyllableData(j,4)/1000) * Fs) - (PrePostPadding * Fs));
            EndIndex = round(((IndividualBirds(Bird).AllSyllableData(j,5)/1000) * Fs) + (PrePostPadding * Fs));
            if ((StartIndex < 1) || (EndIndex > length(LogAmplitude)))
                IndividualBirds(Bird).SyllLogAmplitude{end+1} = NaN;
                IndividualBirds(Bird).SyllLogAmplitudeStatus(end+1) = 0;
                IndividualBirds(Bird).RawSyllLogAmplitudeMeanValue(end+1) = NaN;
                continue;
            else
                IndividualBirds(Bird).SyllLogAmplitudeStatus(end+1) = 1;
                IndividualBirds(Bird).SyllLogAmplitude{end+1} = LogAmplitude(StartIndex:EndIndex);
                TempRawSyllLogAmplitudes{end+1} = LogAmplitude(StartIndex:EndIndex);
                % Now get the values for these syllables for the raw segmented
                % boundaries
                StartIndex = round(((IndividualBirds(Bird).AllSyllableData(j,4)/1000) * Fs));
                EndIndex = round(((IndividualBirds(Bird).AllSyllableData(j,5)/1000) * Fs));
                IndividualBirds(Bird).RawSyllLogAmplitudeMeanValue(end+1) = mean(LogAmplitude(StartIndex:EndIndex));
            end
        end
    end
    
    IndividualBirds(Bird).AdjustedSyllLogAmplitudeMeanValue = ones(size(IndividualBirds(Bird).RawSyllLogAmplitudeMeanValue))*NaN;
    
    % Now to find the common pre and post thresholds and then adjust
    % boundaries and calculate mean amplitude values accordingly
    PrePostPaddingIndex = round(PrePostPadding * Fs);
    for Sylls = 1:length(IndividualBirds(Bird).SortedBirdParameters(1).MotifLabels),
%         figure; hold on;
        PreMinValue{Sylls} = [];
        PostMinValue{Sylls} = [];
        for i = 1:length(IndividualBirds(Bird).SyllLogAmplitude),
            if ((IndividualBirds(Bird).SyllLogAmplitudeStatus(i) == 1) && (char(IndividualBirds(Bird).AllSyllableData(i,1)) == IndividualBirds(Bird).SortedBirdParameters(1).MotifLabels(Sylls))) 
                PreMinValue{Sylls}(i) = min(IndividualBirds(Bird).SyllLogAmplitude{i}(1:PrePostPaddingIndex));
                PostMinValue{Sylls}(i) = min(IndividualBirds(Bird).SyllLogAmplitude{i}(length(IndividualBirds(Bird).SyllLogAmplitude{i}) - PrePostPaddingIndex:end));
            else
                PreMinValue{Sylls}(i) = NaN;
                PostMinValue{Sylls}(i) = NaN;
            end
        end

        PreCommonThreshold{Sylls} = prctile(PreMinValue{Sylls}, ThresholdPercentile);
        PostCommonThreshold{Sylls} = prctile(PostMinValue{Sylls}, ThresholdPercentile);
        
        
        % Now adjust boundaries according to new threshold and calculate mean
        % values

        for i = 1:length(IndividualBirds(Bird).SyllLogAmplitude),
            if ((IndividualBirds(Bird).SyllLogAmplitudeStatus(i) == 1) && (char(IndividualBirds(Bird).AllSyllableData(i,1)) == IndividualBirds(Bird).SortedBirdParameters(1).MotifLabels(Sylls))) 
                if (IndividualBirds(Bird).SyllLogAmplitude{i}(PrePostPaddingIndex + 1) > PreCommonThreshold{Sylls})
                   [Val, NewStartIndex] = find(IndividualBirds(Bird).SyllLogAmplitude{i}(1:PrePostPaddingIndex) <= PreCommonThreshold{Sylls}, 1, 'last');
                else
                   [Val, NewStartIndex] = find(IndividualBirds(Bird).SyllLogAmplitude{i}(PrePostPaddingIndex+1:end) >= PreCommonThreshold{Sylls}, 1, 'first');
                   if (~isempty(NewStartIndex))
                       NewStartIndex = NewStartIndex + PrePostPaddingIndex;
                   end
                end
                if isempty(NewStartIndex)
                   NewStartIndex = PrePostPaddingIndex + 1;
                end

                if (IndividualBirds(Bird).SyllLogAmplitude{i}(end - PrePostPaddingIndex) > PostCommonThreshold{Sylls})
                   [Val, NewEndIndex] = find(IndividualBirds(Bird).SyllLogAmplitude{i}(end - PrePostPaddingIndex + 1:end) <= PostCommonThreshold{Sylls}, 1, 'first');
                   if (~isempty(NewEndIndex))
                       NewEndIndex = NewEndIndex + length(IndividualBirds(Bird).SyllLogAmplitude{i}) - PrePostPaddingIndex;
                   end
                else
                   [Val, NewEndIndex] = find(IndividualBirds(Bird).SyllLogAmplitude{i}(1:end - PrePostPaddingIndex) >= PostCommonThreshold{Sylls}, 1, 'last');
                end
                if isempty(NewEndIndex)
                   NewEndIndex = length(IndividualBirds(Bird).SyllLogAmplitude{i}) - PrePostPaddingIndex;
                end
%                 plot(linspace(1, length(IndividualBirds(Bird).SyllLogAmplitude{i}), length(IndividualBirds(Bird).SyllLogAmplitude{i})/50), spline(1:1:length(IndividualBirds(Bird).SyllLogAmplitude{i}), IndividualBirds(Bird).SyllLogAmplitude{i}, linspace(1, length(IndividualBirds(Bird).SyllLogAmplitude{i}), length(IndividualBirds(Bird).SyllLogAmplitude{i})/50)), 'b', 'LineWidth', 0.1);
%                 axis tight;
%                 Temp = axis;
%                 plot(ones(1,2)*PrePostPaddingIndex, Temp(3:4), 'k--');
%                 plot(ones(1,2)*(Temp(2) - PrePostPaddingIndex), Temp(3:4), 'k--');
                IndividualBirds(Bird).NewStartIndex(i) = NewStartIndex;
                IndividualBirds(Bird).NewEndIndex(i) = NewEndIndex;
                IndividualBirds(Bird).AdjustedSyllLogAmplitudeMeanValue(i) = mean(IndividualBirds(Bird).SyllLogAmplitude{i}(NewStartIndex:NewEndIndex));
                TempAdjustedSyllLogAmplitudes{i} = IndividualBirds(Bird).SyllLogAmplitude{i}(NewStartIndex:NewEndIndex);
            end
        end
%         axis tight;
%         Temp = axis;
%         plot([Temp(1) PrePostPaddingIndex], ones(1,2)*PreCommonThreshold{Sylls}, 'k--');
%         plot([(Temp(2) - PrePostPaddingIndex + 1) Temp(2)], ones(1,2)*PostCommonThreshold{Sylls}, 'k--');

        for i = 1:length(IndividualBirds(Bird).SyllLogAmplitude),
            if ((IndividualBirds(Bird).SyllLogAmplitudeStatus(i) == 1) && (char(IndividualBirds(Bird).AllSyllableData(i,1)) == IndividualBirds(Bird).SortedBirdParameters(1).MotifLabels(Sylls))) 
%                 axis tight;
%                 Temp = axis;
%                 plot(IndividualBirds(Bird).NewStartIndex(i), IndividualBirds(Bird).SyllLogAmplitude{i}(IndividualBirds(Bird).NewStartIndex(i)), 'rs', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
%                 plot(IndividualBirds(Bird).NewEndIndex(i), IndividualBirds(Bird).SyllLogAmplitude{i}(IndividualBirds(Bird).NewEndIndex(i)), 'rs', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
            end
        end
    end
    if (~iscell(BirdNames))
        save(['/home/raghav/StudentRelated/Harini/', BirdNames, 'RawSyllLogAmpTraces.mat'], 'TempRawSyllLogAmplitudes', 'TempAdjustedSyllLogAmplitudes');
    else
        save(['/home/raghav/StudentRelated/Harini/', BirdNames{Bird}, 'RawSyllLogAmpTraces.mat'], 'TempRawSyllLogAmplitudes', 'TempAdjustedSyllLogAmplitudes');
    end

    IndividualBirds(Bird).SyllLogAmplitude = [];
end

RawSyllLogAmplitudeMeanValue = IndividualBirds(1).RawSyllLogAmplitudeMeanValue;
AdjustedSyllLogAmplitudeMeanValue = IndividualBirds(1).AdjustedSyllLogAmplitudeMeanValue;

fprintf('\n');
disp('Done with adjusting log amplitude calculations');