function [] = ASSL_GenerateSummaryReport(handles)

% This function is to generate a summary of the labelling. I basically want
% the following in the summary
% a) spectrogram with examples of all syllable labels
% b) # of occurences of each of the syllables
% c) cluster variance as measured by mean distance of individual points
% from cluster centroid - Mahalanobis distance (only if there are more than
% 5 points)
% d) spectrograms of some of the outliers, chosen based on the criteria
% that Mahalanobis distance of that point is > 75th percentile + 3*IQR
% don't use ff, log amplitude and rms amplitude for this
% e) also show the most common sequences

figure;

Index = 1;
for i = 1:length(handles.ASSL.SyllLabels),
    for j = 1:length(handles.ASSL.SyllLabels{i}),
        AllSylls_Labels(Index) = handles.ASSL.SyllLabels{i}(j);
        AllSylls_FileNo(Index) = i;
        AllSylls_OnsetTime(Index) = handles.ASSL.SyllOnsets{i}(j);
        AllSylls_OffsetTime(Index) = handles.ASSL.SyllOffsets{i}(j);
        for k = 1:length(handles.ASSL.ToBeUsedFeatures),
            if (isempty(strfind(handles.ASSL.ToBeUsedFeatures{k}, 'FundamentalFrequency')))
                AllSylls_SAPFeatures(Index,k) = eval(['handles.ASSL.', handles.ASSL.ToBeUsedFeatures{k}, '{', num2str(i), '}', '(', num2str(j), ')']);
            end
        end
        Index = Index + 1;
    end
end

UniqueSylls = unique(AllSylls_Labels);
for i = 1:length(UniqueSylls),
    SyllCount(i) = length(find(AllSylls_Labels == UniqueSylls(i)));
end

FeaturesNotToBeUsed = {'LogAmplitude' 'RMSAmplitude' 'FundamentalFrequency'};
for i = 1:length(FeaturesNotToBeUsed),
    FeaturesNotToBeUsedIndices(i) = find(strcmp(handles.ASSL.ToBeUsedFeatures, FeaturesNotToBeUsed{i}));
end
FeaturesToBeUsedIndices = setdiff(1:1:length(handles.ASSL.ToBeUsedFeatures), FeaturesNotToBeUsedIndices);

% First a spectrogram for all Capital letter sylls
p = panel();
p.pack(ceil(length(find(SyllCount >= 5))/3), 3);

PlotIndex = 0;
CapitalSylls = find(double(UniqueSylls) < 97);
TimeOffset = 0;
SpectrogramData = [];
% TextStringLocOffset = 0.1;

for i = 1:length(UniqueSylls),
    if (SyllCount(i) >= 5)
        PlotIndex = PlotIndex + 1;
        RelevantSylls = find(AllSylls_Labels == UniqueSylls(i));
        % Distances = pdist2(AllSylls_SAPFeatures(RelevantSylls,FeaturesToBeUsedIndices), median(AllSylls_SAPFeatures(RelevantSylls,FeaturesToBeUsedIndices)), 'mahalanobis', cov(AllSylls_SAPFeatures(RelevantSylls,FeaturesToBeUsedIndices)));
        % Now find min distance and use that syllable for the spectrogram
        %[MinVal, MinValIndex] = min(Distances);
        %MedianSyllIndex = RelevantSylls(MinValIndex);
        SampleSylls = randperm(length(RelevantSylls));
        SampleSylls = SampleSylls(1:5);
        SampleSylls = RelevantSylls(SampleSylls);
        Index = 1;
        SpectrogramData = [];
        for j = 1:length(SampleSylls),
            MedianSyllIndex = SampleSylls(j);
            [Data, Fs] = GetData(handles.ASSL.FileDir{AllSylls_FileNo(MedianSyllIndex)}, handles.ASSL.FileName{AllSylls_FileNo(MedianSyllIndex)}, handles.ASSL.FileType, 0);
            Data = Data(:);
            if (Index == 1)
                SpectrogramData = [SpectrogramData; zeros(round(Fs*0.1), 1)];
            end
            OnsetTimeIndex = round((AllSylls_OnsetTime(MedianSyllIndex) - 0.01) * Fs/1000);
            if (OnsetTimeIndex < 1)
                OnsetTimeIndex = 1;
            end
        
            OffsetTimeIndex = round((AllSylls_OffsetTime(MedianSyllIndex) + 0.01) * Fs/1000);
            if (OffsetTimeIndex > length(Data))
                OffsetTimeIndex = length(Data);
            end

            SpectrogramData = [SpectrogramData; Data(OnsetTimeIndex:OffsetTimeIndex)];
            SpectrogramData = [SpectrogramData; zeros(round(Fs*0.1), 1)];
    %         TextString(Index) = UniqueSylls(CapitalSylls(i));
    %         TextStringLoc(Index) = TextStringLocOffset;
    %         TextStringLocOffset = TextStringLocOffset + length(Data(OnsetTimeIndex:OffsetTimeIndex))/Fs + 0.1;
            Index = Index + 1;
        end
        p(ceil(PlotIndex/3), PlotIndex - 3*(ceil(PlotIndex/3) - 1)).select();
        if (~isempty(SpectrogramData))
            PlotSpectrogramInAxis_SongVar(SpectrogramData, (1:1:length(SpectrogramData))/Fs, Fs, gca);
            set(gca, 'XColor', 'w');
            set(gca, 'YColor', 'w');
            title(UniqueSylls(i));
        end
    end
end
% if (~isempty(SpectrogramData))
%     PlotSpectrogramInAxis_SongVar(SpectrogramData, (1:1:length(SpectrogramData))/Fs, Fs, gca);
%     for i = 1:length(TextStringLoc),
%         text(TextStringLoc(i), 7000, TextString(i), 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'b');
%     end
%     set(gca, 'YColor', 'w');
%     title('Upper case labels');
%     CapitalLabelAxisLimits = axis;
% else
%     CaptitalLabelAxisLimits = [0 0 0 0];
% end

% p(1,2).select();
% LowerCaseSylls = find(double(UniqueSylls) >= 97);
% SpectrogramData = [];
% TextStringLocOffset = 0.1;
% Index = 1;
% clear TextString TextStringLoc;
% 
% 
% for i = 1:length(LowerCaseSylls),
%     if (SyllCount(LowerCaseSylls(i)) >= 5)
%         RelevantSylls = find(AllSylls_Labels == UniqueSylls(LowerCaseSylls(i)));
%         Distances = pdist2(AllSylls_SAPFeatures(RelevantSylls,FeaturesToBeUsedIndices), median(AllSylls_SAPFeatures(RelevantSylls,FeaturesToBeUsedIndices)), 'mahalanobis', cov(AllSylls_SAPFeatures(RelevantSylls,FeaturesToBeUsedIndices)));
%         % Now find min distance and use that syllable for the spectrogram
%         [SortedVals, SortedIndices] = sort(Distances);
%         RelevantSylls = RelevantSylls(SortedIndices);
%         % [MinVal, MinValIndex] = min(Distances);
%         MinValIndex = round(0.015*length(Distances));
%         MedianSyllIndex = RelevantSylls(MinValIndex);
%         [Data, Fs] = GetData(handles.ASSL.FileDir{AllSylls_FileNo(MedianSyllIndex)}, handles.ASSL.FileName{AllSylls_FileNo(MedianSyllIndex)}, handles.ASSL.FileType, 0);
%         Data = Data(:);
%         if (Index == 1)
%             SpectrogramData = [SpectrogramData; zeros(round(Fs*0.1), 1)];
%         end
%         OnsetTimeIndex = round((AllSylls_OnsetTime(MedianSyllIndex) - 0.01) * Fs/1000);
%         if (OnsetTimeIndex < 1)
%             OnsetTimeIndex = 1;
%         end
%         
%         OffsetTimeIndex = round((AllSylls_OffsetTime(MedianSyllIndex) + 0.01) * Fs/1000);
%         if (OffsetTimeIndex > length(Data))
%             OffsetTimeIndex = length(Data);
%         end
%         
%         SpectrogramData = [SpectrogramData; Data(OnsetTimeIndex:OffsetTimeIndex)];
%         SpectrogramData = [SpectrogramData; zeros(round(Fs*0.1), 1)];
%         TextString(Index) = UniqueSylls(LowerCaseSylls(i));
%         TextStringLoc(Index) = TextStringLocOffset;
%         TextStringLocOffset = TextStringLocOffset + length(Data(OnsetTimeIndex:OffsetTimeIndex))/Fs + 0.1;
%         Index = Index + 1;
%     end
% end
% PlotSpectrogramInAxis_SongVar(SpectrogramData, (1:1:length(SpectrogramData))/Fs, Fs, gca);
% for i = 1:length(TextStringLoc),
%     text(TextStringLoc(i), 7000, TextString(i), 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'b');
% end
% set(gca, 'YColor', 'w');
% title('Lower case labels');
% LowerCaseLabelAxisLimits = axis;
% % EqualScaleAxisLimits = [0 max(LowerCaseLabelAxisLimits(2), CapitalLabelAxisLimits(2)) 300 8000];
% % axis(EqualScaleAxisLimits);
% % 
% % p(1,1).select();
% % axis(EqualScaleAxisLimits);
% 
% % Now find most common sequence of syllables
% % I'm going to find this by getting transition probabilities
% for i = 1:length(UniqueSylls),
%     for j = 1:length(UniqueSylls),
%         FirstSylls = find(AllSylls_Labels(1:end-1) == UniqueSylls(i));
%         NextSylls = find(AllSylls_Labels(FirstSylls + 1) == UniqueSylls(j));
%         TransProb(i,j) = length(NextSylls)/length(FirstSylls);
%     end
% end
% 
% % Now find all Transition probabilities > 0.7 with syllable count > 10
% Index = 1;
% for i = 1:length(UniqueSylls),
%     if (SyllCount(i) >= 10)
%         [MaxTransProb, MaxTransProbIndex] = max(TransProb(i,:));
%         if (MaxTransProb >= 0.3)
%             Bigram(Index,:) = [UniqueSylls(i), UniqueSylls(MaxTransProbIndex)];
%             Index = Index + 1;
%         end
%     end
% end
% 
% % Now find longest bigrams
% Index = 1;
% Flag = 1;
% for i = 1:length(Bigram),
%     CurrentBigramIndex = i;
%     CurrentBigram = Bigram(i,:);
%     Flag = 1;
%     Iterations = 0;
%     while ((Flag == 1) && (Iterations < 8))
%         NextBigram = find(Bigram(:,1) == CurrentBigram(end));
%         if (~isempty(NextBigram) && (NextBigram ~= CurrentBigramIndex))
%             CurrentBigram = [CurrentBigram Bigram(NextBigram,2)];
%             CurrentBigramIndex = NextBigram;
%         else
%             Flag = 0;
%         end
%         Iterations = Iterations + 1;
%     end
%     LongestSequence{i} = CurrentBigram;
% end
% MaxSequenceLen = max(cellfun(@length, LongestSequence));
% 
% MaxSequenceIndices = find(cellfun(@length, LongestSequence) == MaxSequenceLen);
% ValidSeqIndex = 1;
% for i = MaxSequenceIndices(:)',
%     Sequences = strfind(AllSylls_Labels, LongestSequence{i});
%     for j = 1:length(Sequences),
%         % Check to see if all members of sequence are in same file
%         SeqOnset = Sequences(j);
%         SeqOffset = Sequences(j) + MaxSequenceLen - 1;
%         SeqOnsetFile = AllSylls_FileNo(SeqOnset);
%         if (SeqOnsetFile ~= mean(AllSylls_FileNo(SeqOnset:SeqOffset)))
%             continue;
%         end
%         
%         % Check to see that the gaps between syllables is not >= 150ms
%         Gaps = AllSylls_OnsetTime(SeqOnset+1:SeqOffset) - AllSylls_OffsetTime(SeqOnset:SeqOffset-1);
%         if (~isempty(find(Gaps >= 150)))
%             continue;
%         end
%         
%         ValidSeqsFileNo(ValidSeqIndex) = AllSylls_FileNo(SeqOnset);
%         ValidSeqsSequence{ValidSeqIndex} = LongestSequence{i};
%         ValidSeqsSequenceOnsetTime(ValidSeqIndex) = AllSylls_OnsetTime(SeqOnset);
%         ValidSeqsSequenceOffsetTime(ValidSeqIndex) = AllSylls_OffsetTime(SeqOffset);
%         ValidSeqIndex = ValidSeqIndex + 1;
%         break;
%     end
% end
% 
% [UniqueValidSeqsFileNo, UniqueFileIndices] = unique(ValidSeqsFileNo);

% p(2).pack('h', {1/2 1/2});
% p(2,1).pack(length(UniqueFileIndices), 1);
% for i = 1:length(UniqueFileIndices),
%     p(2,1,i,1).select();
%     PlotSpectrogramInAxis(handles.ASSL.FileDir{ValidSeqsFileNo(UniqueFileIndices(i))}, handles.ASSL.FileName{ValidSeqsFileNo(UniqueFileIndices(i))}, handles.ASSL.FileType, gca, [ValidSeqsSequenceOnsetTime(UniqueFileIndices(i))/1000 ValidSeqsSequenceOffsetTime(UniqueFileIndices(i))/1000]);
%     set(gca, 'YColor', 'w');
% end

set(gcf, 'Position', [680 289 1100 725]);
set(gcf, 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');

p.margintop = 10;
p.fontsize = 14;

print(fullfile(handles.ASSL.StartDir, 'ASSLNoteFiles', [handles.ASSL.FileListName, '.ASSL.Summary.png']), '-dpng', '-r300');