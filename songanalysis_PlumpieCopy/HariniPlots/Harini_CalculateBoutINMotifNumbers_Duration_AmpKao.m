function [Bout_Stats] = Harini_CalculateBoutINMotifNumbers_Duration_AmpKao(IndividualBirds, Bouts, MotifLabels, INLabels, CommonMotifs)

Bout_Stats.TotalINNumber = [];
Bout_Stats.TotalINNumber_500ms = [];
Bout_Stats.INNumber_500ms_Time = [];
Bout_Stats.TimeToMotif = [];
Bout_Stats.CompleteMotifNumber = [];
Bout_Stats.PartialMotifNumber = [];
Bout_Stats.AllMotifDuration = [];
Bout_Stats.FirstMotifDuration = [];
Bout_Stats.BoutLength = [];
Bout_Stats.RecordingDay = [];
Bout_Stats.RecordingDayIndices = [];

FFSylls = intersect(MotifLabels, IndividualBirds.FF_SyllLabels);
for j = 1:length(FFSylls),
    Bout_Stats.FF_SyllLabel(j) = FFSylls(j);
    Bout_Stats.FF{j} = [];
    Bout_Stats.FF_SyllDuration{j} = [];
    Bout_Stats.FF_SyllOnsets{j} = [];
    Bout_Stats.FF_SyllOffsets{j} = [];
end

LogAmpSylls = intersect(MotifLabels, IndividualBirds.UniqueSyllLabels);
for j = 1:length(LogAmpSylls),
    Bout_Stats.LogAmplitude_SyllLabel(j) = LogAmpSylls(j);
    Bout_Stats.LogAmplitude{j} = [];
    Bout_Stats.SyllDuration{j} = [];
    Bout_Stats.LogAmp_SyllOnsets{j} = [];
    Bout_Stats.LogAmp_SyllOffsets{j} = [];
    Bout_Stats.LogAmp_BoutIndex{j} = [];
end

Index = 0;
BoutIndexColumnIndex = find(cellfun(@length,strfind(IndividualBirds.BoutStatisticsColumnNames, 'BoutIndex')));

for k = Bouts(:)',
    Index = Index + 1;
    BoutLabels = char(IndividualBirds.AllSyllableData(IndividualBirds.Bouts(k,1):IndividualBirds.Bouts(k,2), 1));
    BoutOnsets = IndividualBirds.AllSyllableData(IndividualBirds.Bouts(k,1):IndividualBirds.Bouts(k,2), 4);
    BoutOffsets = IndividualBirds.AllSyllableData(IndividualBirds.Bouts(k,1):IndividualBirds.Bouts(k,2), 5);
    
    FreqColIndex = find(cellfun(@length, strfind(IndividualBirds.SortedBirdParameters(1).SAPFeat_FieldNames, 'MeanFundamentalFrequency')));
    BoutSAPFF = IndividualBirds.AllSyllableFeatValues(IndividualBirds.Bouts(k,1):IndividualBirds.Bouts(k,2),FreqColIndex);
    BoutSAPFF = BoutSAPFF(:);
    
    BoutSAPLogAmp = IndividualBirds.AllSyllableLogAmplitudeKao(IndividualBirds.Bouts(k,1):IndividualBirds.Bouts(k,2));
    BoutSAPLogAmp = BoutSAPLogAmp(:);
    
    DurationColIndex = find(cellfun(@length, strfind(IndividualBirds.SortedBirdParameters(1).SAPFeat_FieldNames, 'Duration')));
    BoutSAPDur = IndividualBirds.AllSyllableFeatValues(IndividualBirds.Bouts(k,1):IndividualBirds.Bouts(k,2),DurationColIndex);
    BoutSAPDur = BoutSAPDur(:);
    
    Bout_Stats.BoutLabels{Index} = BoutLabels(:)';
    Bout_Stats.BoutIndex(Index) = k;
    Bout_Stats.BoutStatisticsBoutIndex(Index) = find(IndividualBirds.BoutStatistics(:,BoutIndexColumnIndex) == k);
    Bout_Stats.RecordingDayIndices(Index) = IndividualBirds.AllRecordingDayIndices(IndividualBirds.Bouts(k,1));
    Bout_Stats.RecordingDay{Index} = IndividualBirds.RecordingDays{Bout_Stats.RecordingDayIndices(Index)};
    
    FirstMotifSyll = NaN;
    for BoutSyll = 1:length(BoutLabels),
        if (~isempty(find(MotifLabels == BoutLabels(BoutSyll))))
            FirstMotifSyll = BoutSyll;
            break;
        end
    end

    if (~isnan(FirstMotifSyll))
        Bout_Stats.TimeToMotif(Index) = BoutOnsets(FirstMotifSyll) - BoutOnsets(1);
        PreMotifSyllIdentity = zeros(1, (FirstMotifSyll - 1));
        for PreMotifSylls = 1:(FirstMotifSyll-1),
            if (~isempty(find(INLabels == BoutLabels(PreMotifSylls))))
                PreMotifSyllIdentity(PreMotifSylls) = 1;
            end
        end
        % Now to find the last set of INs before the first Motif Syll
        % Find last non-IN syllable
        LastNonINSyll = find(PreMotifSyllIdentity == 0, 1, 'last');
        if (isempty(LastNonINSyll))
            INs = 1:1:length(PreMotifSyllIdentity);
        else
            INs = (LastNonINSyll+1):1:length(PreMotifSyllIdentity);
        end
        
        if (~isempty(INs))
            GapsBetweenINs = BoutOnsets(INs(2:end)) - BoutOffsets(INs(1:end-1));
            LongGaps = find(GapsBetweenINs > 500);
            if (isempty(LongGaps))
                INs_With500ms = INs;
            else
                INs_With500ms = INs(LongGaps(end)+1:end);
            end
        else
            INs_With500ms = [];
        end
    else
        Bout_Stats.TimeToMotif(Index) = NaN;
        INs = NaN;
        INs_With500ms = NaN;
    end
    
    CompleteMotifs = strfind(BoutLabels(:)', CommonMotifs{1});
    PartialMotifs = strfind(BoutLabels(:)', CommonMotifs{1}(1:round(length(CommonMotifs{1})/2)));
    if (~isempty(CompleteMotifs))
        MotifDurations = BoutOffsets(CompleteMotifs + length(CommonMotifs{1}) - 1) - BoutOnsets(CompleteMotifs);
        FirstMotifDuration = MotifDurations(1);
    end
    
    if (~isnan(INs))
        Bout_Stats.TotalINNumber(Index) = length(INs);
        Bout_Stats.TotalINNumber_500ms(Index) = length(INs_With500ms);
        if (~isempty(INs))
            Bout_Stats.INNumber_500ms_Time(Index) = BoutOnsets(INs(end)+1) - BoutOnsets(INs(1));
        else
            Bout_Stats.INNumber_500ms_Time(Index) = 0;
        end
    else
        Bout_Stats.TotalINNumber(Index) = NaN;
        Bout_Stats.TotalINNumber_500ms(Index) = NaN;
        Bout_Stats.INNumber_500ms_Time(Index) = NaN;
    end
    
    if (~isempty(CompleteMotifs))
        Bout_Stats.CompleteMotifNumber(Index) = length(CompleteMotifs);
        Bout_Stats.PartialMotifNumber(Index) = length(PartialMotifs);
        Bout_Stats.AllMotifDuration = [Bout_Stats.AllMotifDuration; MotifDurations(:)];
        Bout_Stats.FirstMotifDuration(Index) = FirstMotifDuration;
    else
        Bout_Stats.CompleteMotifNumber(Index) = NaN;
        Bout_Stats.PartialMotifNumber(Index) = NaN;
    end
    
    Bout_Stats.BoutLength(Index) = BoutOffsets(end) - BoutOnsets(1);
    
    % Now for FF
    for j = 1:length(FFSylls),
        MatchingSylls = find(BoutLabels == FFSylls(j));
        if (~isempty(MatchingSylls))
            Bout_Stats.FF{j} = [Bout_Stats.FF{j}; BoutSAPFF(MatchingSylls(:))];
            Bout_Stats.FF_SyllDuration{j} = [Bout_Stats.FF_SyllDuration{j}; BoutSAPDur(MatchingSylls(:))];
            Bout_Stats.FF_SyllOnsets{j} = [Bout_Stats.FF_SyllOnsets{j}; BoutOnsets(MatchingSylls(:))];
            Bout_Stats.FF_SyllOffsets{j} = [Bout_Stats.FF_SyllOffsets{j}; BoutOffsets(MatchingSylls(:))];
        end
    end
    
    % Now for CV of FF for each individual bout only if complete motifs is
    % greater than 2 indicating that all syllables are present in atleast
    % two copies in that bout to be able to calculate CV
    
    if (length(CompleteMotifs) >= 3)
        clear TempCV;
        if (~isempty(FFSylls))
            for j = 1:length(FFSylls),
                MatchingSylls = find(BoutLabels == FFSylls(j));
                TempCV(j) = std(BoutSAPFF(MatchingSylls))/mean(BoutSAPFF(MatchingSylls));
            end
            Bout_Stats.Bout_FFCV(Index) = mean(TempCV);
        else
            Bout_Stats.Bout_FFCV(Index) = NaN;
        end
    else
        Bout_Stats.Bout_FFCV(Index) = NaN;
    end
    
    % Now for Log amplitude
    for j = 1:length(LogAmpSylls),
        MatchingSylls = find(BoutLabels == LogAmpSylls(j));
        if (~isempty(MatchingSylls))
            % Take all syllables
            Bout_Stats.LogAmplitude{j} = [Bout_Stats.LogAmplitude{j}; BoutSAPLogAmp(MatchingSylls(:))];
            Bout_Stats.SyllDuration{j} = [Bout_Stats.SyllDuration{j}; BoutSAPDur(MatchingSylls(:))];
            Bout_Stats.LogAmp_SyllOnsets{j} = [Bout_Stats.LogAmp_SyllOnsets{j}; BoutOnsets(MatchingSylls(:))];
            Bout_Stats.LogAmp_SyllOffsets{j} = [Bout_Stats.LogAmp_SyllOffsets{j}; BoutOffsets(MatchingSylls(:))];
            Bout_Stats.LogAmp_BoutIndex{j} = [Bout_Stats.LogAmp_BoutIndex{j}; ones(length(MatchingSylls), 1)*k];
            % Take only first syllable of the bout
%             Bout_Stats.LogAmplitude{j} = [Bout_Stats.LogAmplitude{j}; BoutSAPLogAmp(MatchingSylls(1))];
%             Bout_Stats.SyllDuration{j} = [Bout_Stats.SyllDuration{j}; BoutSAPDur(MatchingSylls(1))];
        end
    end
end
disp('Done');
    