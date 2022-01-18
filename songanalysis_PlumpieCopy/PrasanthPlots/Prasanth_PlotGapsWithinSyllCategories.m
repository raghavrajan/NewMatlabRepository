function [OutlierEliminatedMeanDursGaps] = Prasanth_PlotGapsWithinSyllCategories(BirdParameters, MinTrialNo, varargin)

% ================= Inter-syllable interval analysis ======================

% Here, I will pull out all the valid bouts and then pick syllables of the
% same type and calculate the mean intervals between these syllables - they
% have to be within a single bout.
rng('default');

Colours = 'rgbkcm';
Symbols = 'o+d^s<>.';

for i = 1:length(BirdParameters),
   % Find all valid song boutsor 
    % Valid bout is one with >= 2s silence before and after the bout. It
    % can be a song bout or a non-song bout - for this analysis it doesn't
    % matter
    
    ValidBouts = find((BirdParameters(i).Bouts(:,8) > 0) & (BirdParameters(i).Bouts(:,9) > 1));
    disp([BirdParameters(i).BirdName, ': ', num2str(length(ValidBouts)), ' valid song bouts']); 
    
    MotifDurs{i} = [];
    MotifGaps{i} = [];
    
    INDurs{i} = [];
    INGaps{i} = [];
    
    ShortCallDurs{i} = [];
    ShortCallGaps{i} = [];
    
    LongCallDurs{i} = [];
    LongCallGaps{i} = [];
    
    First3Durs{i} = [];
    First3Gaps{i} = [];
    
    RestDurs{i} = [];
    RestGaps{i} = [];
    
    AllGaps{i} = [];
    AllDurs{i} = [];
    
    if (nargin > 2)
        for Reps = 1:10000,
            RandFirst3Durs{Reps} = [];
            RandFirst3Gaps{Reps} = [];
            
            RandRestDurs{Reps} = [];
            RandRestGaps{Reps} = [];
        end
    end
    for j = ValidBouts(:)',
        Indices = find(BirdParameters(i).SyllableListBoutNum == j);
        TempBoutLabels = char(BirdParameters(i).SyllableData(Indices,1));
        BoutOnsets = (BirdParameters(i).SyllableData(Indices,4));
        BoutOffsets = (BirdParameters(i).SyllableData(Indices,5));
        
        % Now to recode this in terms of 0s and 1s for motif sylls, INs,
        % short calls and long calls
        
        INs = zeros(size(TempBoutLabels));
        Motifs = zeros(size(TempBoutLabels));
        ShortCalls = zeros(size(TempBoutLabels));
        LongCalls = zeros(size(TempBoutLabels));
        
        for k = 1:length(TempBoutLabels),
            if (~isempty(find(BirdParameters(i).MotifLabels == TempBoutLabels(k))))
                Motifs(k) = 1;
            else
                if (~isempty(find(BirdParameters(i).INLabels == TempBoutLabels(k))))
                    INs(k) = 1;
                else
                    if (~isempty(find(BirdParameters(i).ShortcallLabels == TempBoutLabels(k))))
                        ShortCalls(k) = 1;
                    else
                        if (~isempty(find(BirdParameters(i).LongcallLabels == TempBoutLabels(k))))
                            LongCalls(k) = 1;
                        end
                    end
                end
            end
        end
        
        % Now in each of these categories, I have to check for continuous
        % syllables of that category
        % First for motifs
        MotifOnsets = find(conv(Motifs, [1 -1]) == 1);
        MotifOffsets = find(conv(Motifs, [1 -1]) == -1) - 1;
        if (~isempty(MotifOnsets))
            MotifLengths = MotifOffsets - MotifOnsets;
            ContinuousMotifSylls = find(MotifLengths > 0);
            for k = ContinuousMotifSylls(:)',
                Gaps = BoutOnsets(MotifOnsets(k)+1:MotifOffsets(k)) - BoutOffsets(MotifOnsets(k):MotifOffsets(k)-1);
                Durs = BoutOffsets(MotifOnsets(k):MotifOffsets(k)) - BoutOnsets(MotifOnsets(k):MotifOffsets(k));
                MotifDurs{i} = [MotifDurs{i}; Durs];
                MotifGaps{i} = [MotifGaps{i}; Gaps];
            end
        end
        
        % Next for INs
        INOnsets = find(conv(INs, [1 -1]) == 1);
        INOffsets = find(conv(INs, [1 -1]) == -1) - 1;
        if (~isempty(INOnsets))
            INLengths = INOffsets - INOnsets;
            ContinuousINSylls = find(INLengths > 0);
            for k = ContinuousINSylls(:)',
                Gaps = BoutOnsets(INOnsets(k)+1:INOffsets(k)) - BoutOffsets(INOnsets(k):INOffsets(k)-1);
                Durs = BoutOffsets(INOnsets(k):INOffsets(k)) - BoutOnsets(INOnsets(k):INOffsets(k));
                INDurs{i} = [INDurs{i}; Durs];
                INGaps{i} = [INGaps{i}; Gaps];
            end
        end
        
        % Next for short calls
        ShortCallOnsets = find(conv(ShortCalls, [1 -1]) == 1);
        ShortCallOffsets = find(conv(ShortCalls, [1 -1]) == -1) - 1;
        if (~isempty(ShortCallOnsets))
            ShortCallLengths = ShortCallOffsets - ShortCallOnsets;
            ContinuousShortCallSylls = find(ShortCallLengths > 0);
            for k = ContinuousShortCallSylls(:)',
                Gaps = BoutOnsets(ShortCallOnsets(k)+1:ShortCallOffsets(k)) - BoutOffsets(ShortCallOnsets(k):ShortCallOffsets(k)-1);
                Durs = BoutOffsets(ShortCallOnsets(k):ShortCallOffsets(k)) - BoutOnsets(ShortCallOnsets(k):ShortCallOffsets(k));
                ShortCallDurs{i} = [ShortCallDurs{i}; Durs];
                ShortCallGaps{i} = [ShortCallGaps{i}; Gaps];
            end
        end
        
        % Next for long calls
        LongCallOnsets = find(conv(LongCalls, [1 -1]) == 1);
        LongCallOffsets = find(conv(LongCalls, [1 -1]) == -1) - 1;
        if (~isempty(LongCallOnsets))
            LongCallLengths = LongCallOffsets - LongCallOnsets;
            ContinuousLongCallSylls = find(LongCallLengths > 0);
            for k = ContinuousLongCallSylls(:)',
                Gaps = BoutOnsets(LongCallOnsets(k)+1:LongCallOffsets(k)) - BoutOffsets(LongCallOnsets(k):LongCallOffsets(k)-1);
                Durs = BoutOffsets(LongCallOnsets(k):LongCallOffsets(k)) - BoutOnsets(LongCallOnsets(k):LongCallOffsets(k));
                LongCallDurs{i} = [LongCallDurs{i}; Durs];
                LongCallGaps{i} = [LongCallGaps{i}; Gaps];
            end
        end
        
        if (length(BoutOnsets) >= 2)
            First3Durs{i} = [First3Durs{i}; (BoutOffsets(1:min(3,length(BoutOffsets))) - BoutOnsets(1:min(3,length(BoutOffsets))))];
            First3Gaps{i} = [First3Gaps{i}; (BoutOnsets(2:min(3,length(BoutOffsets))) - BoutOffsets(1:min(2,length(BoutOffsets))))];
        end
        
        if (length(BoutOnsets) >= 5)
            RestDurs{i} = [RestDurs{i}; (BoutOffsets(4:end) - BoutOnsets(4:end))];
            RestGaps{i} = [RestGaps{i}; (BoutOnsets(5:end) - BoutOffsets(4:end-1))];
        end
        
        if (nargin > 2)
            % Randomise syllable onsets and offsets within the bout and then
            % calculate the First 3 Durs and First 3 Gaps
            for Reps = 1:10000,
                RandDurs = BoutOffsets - BoutOnsets;
                RandDurs = RandDurs(randperm(length(RandDurs)));
                
                RandGaps = BoutOnsets(2:end) - BoutOffsets(1:end-1);
                RandGaps = RandGaps(randperm(length(RandGaps)));

                if (length(RandDurs) >= 2)
                    RandFirst3Durs{Reps} = [RandFirst3Durs{Reps}; (RandDurs(1:min(3,length(RandDurs))))];
                    RandFirst3Gaps{Reps} = [RandFirst3Gaps{Reps}; (RandGaps(1:min(2,length(RandGaps))))];
                end
                if (length(RandDurs) >= 5)
                    RandRestDurs{Reps} = [RandRestDurs{Reps}; (RandDurs(4:end))];
                    RandRestGaps{Reps} = [RandRestGaps{Reps}; (RandGaps(4:end))];
                end
            end
        end
    end
    
    if (nargin > 2)
        for Reps = 1:10000,
            RandMeanDursGaps{Reps}(i,:) = [mean(RandFirst3Durs{Reps}) mean(RandFirst3Gaps{Reps}) mean(RandRestDurs{Reps}) mean(RandRestGaps{Reps})];
        end
    end
end

for i = 1:length(BirdParameters),
    MeanDursGaps(i,:) = [mean(MotifDurs{i}) mean(MotifGaps{i}) mean(INDurs{i}) mean(INGaps{i}) mean(ShortCallDurs{i}) mean(ShortCallGaps{i}) mean(LongCallDurs{i}) mean(LongCallGaps{i}) mean(First3Durs{i}) mean(First3Gaps{i}) mean(RestDurs{i}) mean(RestGaps{i})];
    %MeanDursGaps(i,:) = [mean(MotifDurs{i}) mean(MotifGaps{i}) mean(INDurs{i}) mean(INGaps{i}) mean(ShortCallDurs{i}) mean(ShortCallGaps{i}) mean(LongCallDurs{i}) mean(LongCallGaps{i})];
    MedianDursGaps(i,:) = [median(MotifDurs{i}) median(MotifGaps{i}) median(INDurs{i}) median(INGaps{i}) median(ShortCallDurs{i}) median(ShortCallGaps{i}) median(LongCallDurs{i}) median(LongCallGaps{i})];
end

% This figure plots all the syllable types
MeanFig = figure;
hold on;
%MedianFig = figure;
%hold on;

% For confidence ellipse exclude outliers in each column
OutlierEliminatedMeanDursGaps = MeanDursGaps;
for j = 1:size(MeanDursGaps,2),
    Outliers = find((MeanDursGaps(:,j) < (prctile(MeanDursGaps(:,j), 25) - 3 * iqr(MeanDursGaps(:,j)))) | (MeanDursGaps(:,j) > (prctile(MeanDursGaps(:,j), 75) + 3 * iqr(MeanDursGaps(:,j)))));
    if (~isempty(Outliers))
        OutlierEliminatedMeanDursGaps(Outliers, j) = NaN;
    end
end

for i = 1:size(MeanDursGaps,2)/2,
    figure(MeanFig);
    PlotHandle(i) = plot(MeanDursGaps(:,(i-1)*2 + 1), MeanDursGaps(:,i*2), [Colours(i), 's'], 'MarkerFaceColor', Colours(i));
    
    if (nargin <= 2)
        % Now use only non-NaN rows for ellipse
        [NanRows, NanCols] = find(isnan(OutlierEliminatedMeanDursGaps(:, (i-1)*2 + 1:i*2)));
        NonNanRows = setdiff(1:1:size(OutlierEliminatedMeanDursGaps,1), unique(NanRows));
        PlotConfidenceEllipse(OutlierEliminatedMeanDursGaps(unique(NonNanRows),(i-1)*2 + 1:i*2), Colours(i), 1);
    else
        Inputs = varargin{1};
        [NanRows, NanCols] = find(isnan(Inputs(:, (i-1)*2 + 1:i*2)));
        NonNanRows = setdiff(1:1:size(Inputs,1), unique(NanRows));
        PlotConfidenceEllipse(Inputs(NonNanRows,(i-1)*2 + 1:i*2), Colours(i), 1);
    end
%     figure(MedianFig);
%     plot(MedianDursGaps(:,(i-1)*2 + 1), MedianDursGaps(:,i*2), [Colours(i), 'o'], 'MarkerFaceColor', Colours(i));
%     PlotConfidenceEllipse(MedianDursGaps(:,(i-1)*2 + 1:i*2), Colours(i), 1);
    switch i
        case 1
            PlotLegend{i} = 'Motif syllables';
        case 2
            PlotLegend{i} = 'INs';
        case 3
            PlotLegend{i} = 'Short Calls';
        case 4
            PlotLegend{i} = 'Long Calls';
        case 5
            PlotLegend{i} = 'First 3 syllables';
        case 6
            PlotLegend{i} = 'Rest syllables';
    end
end



legend(PlotHandle, PlotLegend);
xlabel('Mean syllable duration (msec)');
ylabel('Mean inter-syllable interval (msec)');
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 269 600 700]);
axis tight;
Temp = axis;
Temp = [min(0, Temp(1)) 1.02*Temp(2) min(0, Temp(3)) 1.02*Temp(4)];
axis(Temp);

% This figure plots only motif syllables, INs, first 3 sylls and rest 3
% sylls for normal birds and only first 3 sylls and rest 3 sylls for lesion
% birds

MeanFig2 = figure;
hold on;
%MedianFig = figure;
%hold on;

PlotLegend = [];
PlotHandle = [];
for i = [1 2 5 6]
    figure(MeanFig2);
    if (nargin <=2)
        PlotHandle(end+1) = plot(MeanDursGaps(:,(i-1)*2 + 1), MeanDursGaps(:,i*2), [Colours(i), 's'], 'MarkerFaceColor', Colours(i));
    else
        if (i > 2)
            PlotHandle(end+1) = plot(MeanDursGaps(:,(i-1)*2 + 1), MeanDursGaps(:,i*2), [Colours(i), 's'], 'MarkerFaceColor', Colours(i));
        end
    end
    if (nargin <= 2)
        % For confidence ellipse exclude outliers in each column
        OutlierEliminatedMeanDursGaps = MeanDursGaps;
        for j = 1:size(MeanDursGaps,2),
            Outliers = find((MeanDursGaps(:,j) < (prctile(MeanDursGaps(:,j), 25) - 3 * iqr(MeanDursGaps(:,j)))) | (MeanDursGaps(:,j) > (prctile(MeanDursGaps(:,j), 75) + 3 * iqr(MeanDursGaps(:,j)))));
            if (~isempty(Outliers))
                OutlierEliminatedMeanDursGaps(Outliers, j) = NaN;
            end
        end
        % Now use only non-NaN rows for ellipse
        [NanRows, NanCols] = find(isnan(OutlierEliminatedMeanDursGaps(:, (i-1)*2 + 1:i*2)));
        NonNanRows = setdiff(1:1:size(OutlierEliminatedMeanDursGaps,1), unique(NanRows));
        PlotConfidenceEllipse(OutlierEliminatedMeanDursGaps(unique(NonNanRows),(i-1)*2 + 1:i*2), Colours(i), 1);
    else
        if (i < 5)
            continue;
        end
        Inputs = varargin{1};
        [NanRows, NanCols] = find(isnan(Inputs(:, (i-1)*2 + 1:i*2)));
        NonNanRows = setdiff(1:1:size(Inputs,1), unique(NanRows));
        PlotConfidenceEllipse(Inputs(NonNanRows,(i-1)*2 + 1:i*2), Colours(i), 1);
    end
%     figure(MedianFig);
%     plot(MedianDursGaps(:,(i-1)*2 + 1), MedianDursGaps(:,i*2), [Colours(i), 'o'], 'MarkerFaceColor', Colours(i));
%     PlotConfidenceEllipse(MedianDursGaps(:,(i-1)*2 + 1:i*2), Colours(i), 1);
    switch i
        case 1
            if (nargin <= 2)
                PlotLegend{end+1} = 'Motif syllables';
            end
        case 2
            if (nargin <= 2)
                PlotLegend{end+1} = 'INs';
            end
        case 5
            PlotLegend{end+1} = 'First 3 syllables';
        case 6
            PlotLegend{end+1} = 'Rest syllables';
    end
end



legend(PlotHandle, PlotLegend);
xlabel('Mean syllable duration (msec)');
ylabel('Mean inter-syllable interval (msec)');
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 269 600 700]);
axis tight;
Temp = axis;
Temp = [min(0, Temp(1)) 1.02*Temp(2) min(0, Temp(3)) 1.02*Temp(4)];
axis(Temp);

% This figure plots only motif syllables, INs, for normal birds and only 
% first 3 sylls and rest 3 sylls for lesion
% birds. In addition for lesion birds, also include ellipses for the lesion
% bird data

MeanFig3 = figure;
hold on;
%MedianFig = figure;
%hold on;

PlotLegend = [];
PlotHandle = [];
for i = [1 2 5 6]
    figure(MeanFig3);
    if (nargin <=2)
        if (i <= 2)
            PlotHandle(end+1) = plot(MeanDursGaps(:,(i-1)*2 + 1), MeanDursGaps(:,i*2), [Colours(i), 's'], 'MarkerFaceColor', Colours(i));
        end
    else
        if (i > 2)
            PlotHandle(end+1) = plot(MeanDursGaps(:,(i-1)*2 + 1), MeanDursGaps(:,i*2), [Colours(i), 's'], 'MarkerFaceColor', Colours(i));
        end
    end
    if (nargin <= 2)
        if (i > 2)
            continue;
        end
        % For confidence ellipse exclude outliers in each column
        OutlierEliminatedMeanDursGaps = MeanDursGaps;
        for j = 1:size(MeanDursGaps,2),
            Outliers = find((MeanDursGaps(:,j) < (prctile(MeanDursGaps(:,j), 25) - 3 * iqr(MeanDursGaps(:,j)))) | (MeanDursGaps(:,j) > (prctile(MeanDursGaps(:,j), 75) + 3 * iqr(MeanDursGaps(:,j)))));
            if (~isempty(Outliers))
                OutlierEliminatedMeanDursGaps(Outliers, j) = NaN;
            end
        end
        % Now use only non-NaN rows for ellipse
        [NanRows, NanCols] = find(isnan(OutlierEliminatedMeanDursGaps(:, (i-1)*2 + 1:i*2)));
        NonNanRows = setdiff(1:1:size(OutlierEliminatedMeanDursGaps,1), unique(NanRows));
        PlotConfidenceEllipse(OutlierEliminatedMeanDursGaps(unique(NonNanRows),(i-1)*2 + 1:i*2), Colours(i), 1);
    else
        if (i < 5)
            continue;
        end
        Inputs = varargin{1};
        [NanRows, NanCols] = find(isnan(Inputs(:, (i-1)*2 + 1:i*2)));
        NonNanRows = setdiff(1:1:size(Inputs,1), unique(NanRows));
        PlotConfidenceEllipse(Inputs(NonNanRows,(i-1)*2 + 1:i*2), Colours(i), 1);
        
        [NanRows, NanCols] = find(isnan(OutlierEliminatedMeanDursGaps(:, (i-1)*2 + 1:i*2)));
        NonNanRows = setdiff(1:1:size(OutlierEliminatedMeanDursGaps,1), unique(NanRows));
        if (i == 5)
            PlotConfidenceEllipse_Boundary(OutlierEliminatedMeanDursGaps(NonNanRows,(i-1)*2 + 1:i*2), Colours(i), 1);
        else
            PlotConfidenceEllipse_Boundary(OutlierEliminatedMeanDursGaps(NonNanRows,(i-1)*2 + 1:i*2), Colours(i), 1);
        end
    end
%     figure(MedianFig);
%     plot(MedianDursGaps(:,(i-1)*2 + 1), MedianDursGaps(:,i*2), [Colours(i), 'o'], 'MarkerFaceColor', Colours(i));
%     PlotConfidenceEllipse(MedianDursGaps(:,(i-1)*2 + 1:i*2), Colours(i), 1);
    switch i
        case 1
            PlotLegend{end+1} = 'Motif syllables';
        case 2
            PlotLegend{end+1} = 'INs';
        case 5
            PlotLegend{end+1} = 'First 3 syllables';
        case 6
            PlotLegend{end+1} = 'Rest syllables';
    end
end



legend(PlotHandle, PlotLegend);
xlabel('Mean syllable duration (msec)');
ylabel('Mean inter-syllable interval (msec)');
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 269 600 700]);
axis tight;
Temp = axis;
Temp = [min(0, Temp(1)) 1.02*Temp(2) min(0, Temp(3)) 1.02*Temp(4)];
axis(Temp);

% Now calculate mahalanobis distance of all points from the First 3
% position and Rest position ellipses of the normal bird
First3Durs_NormalBirds = Inputs(:,(5-1)*2 + 1:5*2);
RestDurs_NormalBirds = Inputs(:,(6-1)*2 + 1:6*2);

First3Distance_First3 = mean(pdist2(OutlierEliminatedMeanDursGaps(:,9:10), mean(First3Durs_NormalBirds), 'mahalanobis', cov(First3Durs_NormalBirds)));
First3Distance_Rest = mean(pdist2(OutlierEliminatedMeanDursGaps(:,9:10), mean(RestDurs_NormalBirds), 'mahalanobis', cov(RestDurs_NormalBirds)));

RestDistance_First3 = mean(pdist2(OutlierEliminatedMeanDursGaps(:,11:12), mean(First3Durs_NormalBirds), 'mahalanobis', cov(First3Durs_NormalBirds)));
RestDistance_Rest = mean(pdist2(OutlierEliminatedMeanDursGaps(:,11:12), mean(RestDurs_NormalBirds), 'mahalanobis', cov(RestDurs_NormalBirds)));

% Now eliminate outliers for the random values of durs and gaps
for i = 1:10000,
    for j = 1:size(RandMeanDursGaps{i}, 2),
        Outliers = find((RandMeanDursGaps{i}(:,j) < (prctile(RandMeanDursGaps{i}(:,j), 25) - 3*iqr(RandMeanDursGaps{i}(:,j)))) | (RandMeanDursGaps{i}(:,j) < (prctile(RandMeanDursGaps{i}(:,j), 75) + 3*iqr(RandMeanDursGaps{i}(:,j)))));
        if (isempty(Outliers))
            RandMeanDursGaps{i}(Outliers,j) = NaN;
        end
    end
    [NanRows, NanCols] = find(isnan(RandMeanDursGaps{i}));
    RandMeanDursGaps(unique(NanRows),:) = [];
end
    
% Now for the chance
for i = 1:10000,
    RandFirst3Distance_First3(i) = mean(pdist2(RandMeanDursGaps{i}(:,1:2), mean(First3Durs_NormalBirds), 'mahalanobis', cov(First3Durs_NormalBirds)));
    RandFirst3Distance_Rest(i) = mean(pdist2(RandMeanDursGaps{i}(:,1:2), mean(RestDurs_NormalBirds), 'mahalanobis', cov(RestDurs_NormalBirds)));

    RandRestDistance_First3(i) = mean(pdist2(RandMeanDursGaps{i}(:,3:4), mean(First3Durs_NormalBirds), 'mahalanobis', cov(First3Durs_NormalBirds)));
    RandRestDistance_Rest(i) = mean(pdist2(RandMeanDursGaps{i}(:,3:4), mean(RestDurs_NormalBirds), 'mahalanobis', cov(RestDurs_NormalBirds)));
end

figure;
hold on;
plot([1 2 4 5], [First3Distance_First3 First3Distance_Rest RestDistance_First3 RestDistance_Rest], 'ks-');
plot([1 2 4 5], [prctile(RandFirst3Distance_First3, 2.5) prctile(RandFirst3Distance_Rest, 2.5) prctile(RandRestDistance_First3, 2.5) prctile(RandRestDistance_Rest, 2.5)], 'k--'); 
plot([1 2 4 5], [prctile(RandFirst3Distance_First3, 97.5) prctile(RandFirst3Distance_Rest, 97.5) prctile(RandRestDistance_First3, 97.5) prctile(RandRestDistance_Rest, 97.5)], 'k--');
set(gca, 'XTick', [1 2 4 5], 'XTickLabel', [{'From First 3'} {'From Rest'} {'From First 3'} {'From Rest'}]); 
ylabel('Mean Distance from centre of cluster');
disp('Finished');