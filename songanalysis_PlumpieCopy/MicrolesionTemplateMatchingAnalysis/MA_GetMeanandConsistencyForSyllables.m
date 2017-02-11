function [Results] = MA_GetMeanandConsistencyForSyllables(Data, Threshold, BoutLengths)

% First get all template match values from Data
TemplateMatchVals = [];
for j = 1:length(Data),
    if (~isempty(Data{j}))
        TemplateMatchVals = [TemplateMatchVals; Data{j}];
    end
end

% Now get mean, median, std, and sem of match values above threshold
Indices = find(TemplateMatchVals(:,1) >= Threshold);
if (~isempty(Indices))
    Results.MeanMatchValue = [mean(TemplateMatchVals(Indices,1)) median(TemplateMatchVals(Indices,1)) std(TemplateMatchVals(Indices,1)) (std(TemplateMatchVals(Indices,1))/sqrt(length(Indices)))];
else
    Results.MeanMatchValue = [NaN NaN NaN NaN];
end
% Now I have to identify which bouts these match values fall within and
% arrange them according to their bout
if (~isempty(Indices))
    for i = 1:length(Indices),
        for j = 1:size(BoutLengths{TemplateMatchVals(Indices(i),3)},1),
            BoutOnset = BoutLengths{TemplateMatchVals(Indices(i),3)}(j,2)/1000;
            BoutOffset = BoutLengths{TemplateMatchVals(Indices(i),3)}(j,3)/1000;
            if ((TemplateMatchVals(Indices(i),2) >= BoutOnset) && (TemplateMatchVals(Indices(i),2) <= BoutOffset))
                Results.MatchVals(i,:) = [TemplateMatchVals(Indices(i),:) j];
            end
        end
    end
else
    Results.MatchVals = [];
end

% Now get the number of match values above threshold per bout divided by 
% bout length - this is a measure of consistency
Results.IntervalBetweenMatches = [];

if (~isempty(Results.MatchVals))
    Index = 1;
    for i = 1:length(BoutLengths),
        for j = 1:size(BoutLengths{i},1),
            Indices = find((Results.MatchVals(:,3) == i) & (Results.MatchVals(:,4) == j));
            Results.MeanNumPerBout(Index) = length(Indices)/(BoutLengths{i}(j,4)/1000);
            if (length(Indices) > 1)
                Results.IntervalBetweenMatches = [Results.IntervalBetweenMatches; diff(Results.MatchVals(Indices,2))];
            end
            Index = Index + 1;
        end
    end
else
    Results.MeanNumPerBout = [];
end

if (~isempty(Results.MatchVals))
    Results.Consistency = [mean(Results.MeanNumPerBout) median(Results.MeanNumPerBout) std(Results.MeanNumPerBout) std(Results.MeanNumPerBout)/sqrt(length(Results.MeanNumPerBout))];
    if (length(Results.IntervalBetweenMatches) > 1)
        Results.MeanIntervalBetweenMatches = [mean(Results.IntervalBetweenMatches) median(Results.IntervalBetweenMatches) std(Results.IntervalBetweenMatches) std(Results.IntervalBetweenMatches)/length(Results.IntervalBetweenMatches)];    
    else
        Results.MeanIntervalBetweenMatches = [];
    end
else
    Results.Consistency = [0 0 0 0];
    Results.MeanIntervalBetweenMatches = [];
end
