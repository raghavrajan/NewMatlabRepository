function [BaselineSubtractedValues, RawValues] = Harini_CalcBaselineSubtractedLogAmplitude(IndividualBirds, BirdNames)

% First remove outliers based on 3 conditions:% 1. Outlier according to template match value
% 2. Outlier according to mahalanobis distance
% 3. Outlier according to correlation value with mean amplitude waveform

% First according to template match value

% All of the above is only for motif syllables that I want to calculate
% amplitude for

Padding = 0.025;
FFTWindowSize = 0.005; % in seconds
FFTWindowStep = 0.001; % in seconds
PrePostTime = 2000; % in milliseconds

for i = 1:length(IndividualBirds),
    BaselineSubtractedValues(i).MeanFFTLogAmpValues = ones(size(IndividualBirds(i).AllFFTLogAmplitudes))*NaN;
    RawValues(i).MeanFFTLogAmpValues = ones(size(IndividualBirds(i).AllFFTLogAmplitudes))*NaN;
    MotifSylls = IndividualBirds(i).SortedBirdParameters(1).MotifLabels;
    ValidSongBouts = find((IndividualBirds(i).Bouts(:,7) == 1) & (IndividualBirds(i).Bouts(:,8) > 0) & (IndividualBirds(i).Bouts(:,9) > 1));
%     for j = ValidSongBouts(:)',
%         if (IndividualBirds(i).SortedBirdParameters(1).Continuousdata == 0)
%             [DataDir, SongFileName, Extension] = fileparts(IndividualBirds(i).SongFileNames{IndividualBirds(i).Bouts(j,3)});
%             [RawData, Fs] = GetData(DataDir, [SongFileName, Extension], IndividualBirds(i).SortedBirdParameters(1).FileType, 1);
%             [LogAmplitude] = LSINA_CalcFFTLogAmplitude(RawData, Fs, FFTWindowSize, FFTWindowStep);
%             Time = linspace(0, length(RawData)/Fs, length(LogAmplitude));
%             LogAmp_Fs = 1/(Time(2) - Time(1));
% %             figure;
%             StartIndex = round((IndividualBirds(i).Bouts(j,5)-PrePostTime)*LogAmp_Fs/1000);
%             EndIndex = round((IndividualBirds(i).Bouts(j,6)+PrePostTime)*LogAmp_Fs/1000);
% %             plot(Time(StartIndex:EndIndex), LogAmplitude(StartIndex:EndIndex), 'k');
% %             uiwait(gcf);
%         end
%     end
    for j = 1:length(MotifSylls),
        MatchSyllables = find((char(IndividualBirds(i).AllSyllableData(:,1)) == MotifSylls(j)) & (IndividualBirds(i).AllSyllableLogAmpStatus(:) == 1));
        
        % First find outliers based on template match values
%        TemplateMatchOutlierThresholds = [(prctile(IndividualBirds(i).AllSyllableTemplateMatchValues(MatchSyllables), 75) + 3*iqr(IndividualBirds(i).AllSyllableTemplateMatchValues(MatchSyllables))) (prctile(IndividualBirds(i).AllSyllableTemplateMatchValues(MatchSyllables), 25) - 3*iqr(IndividualBirds(i).AllSyllableTemplateMatchValues(MatchSyllables)))];
%        TemplateMatchOutliers = find((IndividualBirds(i).AllSyllableTemplateMatchValues(MatchSyllables) > TemplateMatchOutlierThresholds(1)) | (IndividualBirds(i).AllSyllableTemplateMatchValues(MatchSyllables) < TemplateMatchOutlierThresholds(2)));
%        TemplateMatchOutliers = MatchSyllables(TemplateMatchOutliers);
        
        % Now find outliers based on mahalanobis distance
        FeatureCols = [1 3 4 6 7 8]; % exclude amplitude and amplitude modulation
        
        % Now find outliers based on mahalanobis distance
        Matches = MatchSyllables;
        [NanRows, NanCols] = find(isnan(IndividualBirds(i).AllSyllableFeatValues(Matches,FeatureCols)));
        NanRows = unique(NanRows);
        NanRows = Matches(NanRows);
        NonNanRows = setdiff(Matches, NanRows);
        Distances = pdist2(IndividualBirds(i).AllSyllableFeatValues(NonNanRows,FeatureCols), mean(IndividualBirds(i).AllSyllableFeatValues(NonNanRows,FeatureCols)), 'mahalanobis', cov(IndividualBirds(i).AllSyllableFeatValues(NonNanRows,FeatureCols)));
        
        DistanceOutlierThreshold = [(prctile(Distances, 75) + 3*iqr(Distances)) (prctile(Distances, 25) - 3*iqr(Distances))];

        DistanceOutliers = find((Distances > DistanceOutlierThreshold(1)) | (Distances < DistanceOutlierThreshold(2)));
        DistanceOutliers = NonNanRows(DistanceOutliers);
        
        % Now find outliers based on correlation with average amplitude
        % waveform
        AmpWFLens = cellfun(@length, IndividualBirds(i).AllFFTLogAmplitudes(MatchSyllables));
        clear AllAmpWF;
        for k = 1:length(Matches),
            AllAmpWF(k,1:min(AmpWFLens)) = IndividualBirds(i).AllFFTLogAmplitudes{MatchSyllables(k)}(1:min(AmpWFLens));
        end
        MeanAmpWF = mean(AllAmpWF);
        clear Corr;
        for k = 1:length(MatchSyllables),
            Corr(k) = corr(IndividualBirds(i).AllFFTLogAmplitudes{MatchSyllables(k)}(1:min(AmpWFLens))', MeanAmpWF(:));
        end
        CorrOutliers = find(Corr < (prctile(Corr, 25) - 3*iqr(Corr)));
        CorrOutliers = MatchSyllables(CorrOutliers);
        
%        UniqueOutliers = unique([TemplateMatchOutliers(:); DistanceOutliers(:); CorrOutliers(:)]);
        
        %disp([BirdNames{i}, ': Syll ', MotifSylls(j), ': removed ', num2str(length(UniqueOutliers)), ' out of a total of ', num2str(length(MatchSyllables)), ' (', num2str(100*length(UniqueOutliers)/length(MatchSyllables)), '%)']);
        %MatchSyllables = setdiff(MatchSyllables, UniqueOutliers);
        
        % Initialize mean FFT log amplitude values
        
        clear MinValBeginning MinValEnd;
        % First find the minimum value for the beginning and the ends
        for k = 1:length(MatchSyllables),
            MinValBeginning(k) = min(IndividualBirds(i).AllFFTLogAmplitudes{MatchSyllables(k)}(1:round(1.2*Padding*IndividualBirds(i).SortedBirdParameters(1).Amp_Fs{1}(1))));
            MinValEnd(k) = min(IndividualBirds(i).AllFFTLogAmplitudes{MatchSyllables(k)}((end -  1.2*round(1.2*Padding*IndividualBirds(i).SortedBirdParameters(1).Amp_Fs{1}(1))):end));
        end
        MinValBeginning = max(MinValBeginning);
        MinValEnd = max(MinValEnd);
        
        ThresholdVal = max([MinValBeginning MinValEnd]);
        
        % Now for each amplitude waveform, find the nearest end points that
        % are >= the min vals for beginning and end. Then get the mean
        % between those points as the mean amplitude and subtract the
        % baseline for that file from this value
        % Upsample the waveform 10 fold before checking for > threshold
        for k = 1:length(MatchSyllables)
            Amp_X = linspace(0,length(IndividualBirds(i).AllFFTLogAmplitudes{MatchSyllables(k)})/IndividualBirds(i).SortedBirdParameters(1).Amp_Fs{1}(1), length(IndividualBirds(i).AllFFTLogAmplitudes{MatchSyllables(k)}));
            UpSampledAmplitude = spline(linspace(0,length(IndividualBirds(1).AllFFTLogAmplitudes{i})/IndividualBirds(1).SortedBirdParameters(1).Amp_Fs{1}(1), length(IndividualBirds(1).AllFFTLogAmplitudes{i})), IndividualBirds(1).AllFFTLogAmplitudes{i}, linspace(0,length(IndividualBirds(1).AllFFTLogAmplitudes{i})/IndividualBirds(1).SortedBirdParameters(1).Amp_Fs{1}(1), length(IndividualBirds(1).AllFFTLogAmplitudes{i})*10));
            SyllStartIndex = round(Padding*IndividualBirds(i).SortedBirdParameters(1).Amp_Fs{1}(1));
            SyllEndIndex = length(IndividualBirds(i).AllFFTLogAmplitudes{MatchSyllables(k)}) - SyllStartIndex;
            AboveThreshold = ((IndividualBirds(i).AllFFTLogAmplitudes{MatchSyllables(k)} >= ThresholdVal));
            Starts = find(conv(double(AboveThreshold), [1 -1], 'same') == 1);
            Ends = find(conv(double(AboveThreshold), [1 -1], 'same') == -1);
            
            [SyllStart, SyllStartIndex] = min(abs(Starts - SyllStartIndex));
            SyllStartIndex = Starts(SyllStartIndex);
            
            [SyllEnd, SyllEndIndex] = min(abs(Ends - SyllEndIndex));
            SyllEndIndex = Ends(SyllEndIndex);
            RawValues(i).MeanFFTLogAmpValues(MatchSyllables(k)) = mean(IndividualBirds(i).AllFFTLogAmplitudes{MatchSyllables(k)}(SyllStartIndex:SyllEndIndex));
            BaselineSubtractedValues(i).MeanFFTLogAmpValues(MatchSyllables(k)) = mean(IndividualBirds(i).AllFFTLogAmplitudes{MatchSyllables(k)}(SyllStartIndex:SyllEndIndex)) - IndividualBirds(i).AllSyllableLogAmpBaselineAmpValue(MatchSyllables(k));
        end
    end
end

disp('Done with log amplitude calculations');