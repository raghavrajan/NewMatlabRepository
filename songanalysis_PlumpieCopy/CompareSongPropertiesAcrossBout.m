function [MeanSyllFeats, STDSyllFeats, SEMSyllFeats] = CompareSongPropertiesAcrossBout(INR, MotifSylls)

% I want to compare the acoustic properties of syllables across different
% motifs within a bout


% First get all syllables
UniqueSylls = unique([INR.BoutDetails.labels]);

% Using a 3 column variable for storing data
% First column gives the bout number
% Second column gives the syllable index for the match
% Third column is a number n that means it is the nth occurence of that
% syllable in the bout
% Fourth column is whether there were intro notes in-between or not - 0 if
% there were no intro notes in-between and 1 if there were

for i = 1:length(UniqueSylls),
    SyllFeat{i} = [];
    SyllIndices{i} = [];
    SyllRatioToFirstinBout{i} = [];
    for j = 1:length(INR.BoutDetails),
        
        Matches = find(INR.BoutDetails(j).labels == UniqueSylls(i));
        if (~isempty(Matches))
            SyllIndices{i}(end+1,:) = [j Matches(1) 1 0];
            SyllFeat{i}(end+1,:) = INR.BoutDetails(j).Feats(Matches(1),:);
            SyllRatioToFirstinBout{i}(end+1,:) = SyllFeat{i}(end,:)./INR.BoutDetails(j).Feats(Matches(1),:);

            for k = 2:length(Matches),
                if (isempty(find((INR.BoutDetails(j).labels((Matches(k-1)+1):(Matches(k)-1))) == 'i')))
                    SyllIndices{i}(end+1,:) = [j Matches(k) k 0];
                else
                    SyllIndices{i}(end+1,:) = [j Matches(k) k 1];
                end
                SyllFeat{i}(end+1,:) = INR.BoutDetails(j).Feats(Matches(k),:);
                SyllRatioToFirstinBout{i}(end+1,:) = SyllFeat{i}(end,:)./INR.BoutDetails(j).Feats(Matches(1),:);
            end
        end    
    end
end

% First I need to figure out how many points are there for each motif
% position. Only cases which are greater than a minimum number will be
% considered for analysis. I will do this check individually for each of
% the syllables, as birds need not always sing full motifs.

MinNumber = 5;

for i = 1:length(UniqueSylls),
    if (isempty(find(MotifSylls == UniqueSylls(i))))
        continue;
    end
    
    for k = 1:max(SyllIndices{i}(:,3)),
        NumTrials{i}(k) = length(find(SyllIndices{i}(:,3) == k));
    end
    
    MaxMotifNum = find(NumTrials{i} >= MinNumber, 1, 'last');
    PlotNumRows = 4;
    PlotNumCols = MaxMotifNum;

    figure;
    % Each row in the figure is for one of the features
    % The first set of figures plot the individual features for each position in the bout vs. the
    % features in the first motif
    % The last figure is a bar graph summarising the previous graphs
    
    % First find the maximum bout position
     
    for j = 1:PlotNumRows,
        FirstMotifSylls = find(SyllIndices{i}(:,3) == 1);
        MeanSyllFeats{i}(1) = mean(SyllRatioToFirstinBout{i}(FirstMotifSylls, j));
        STDSyllFeats{i}(1) = std(SyllRatioToFirstinBout{i}(FirstMotifSylls, j));
        SEMSyllFeats{i}(1) = std(SyllRatioToFirstinBout{i}(FirstMotifSylls, j))/sqrt(length(FirstMotifSylls));
        
        for k = 2:MaxMotifNum,
            subplot(PlotNumRows, PlotNumCols, ((j-1)*PlotNumCols) + (k-1));
            hold on;
            Indices = find(SyllIndices{i}(:,3) == k);
            CorrespondingFirstMotifIndices = (Indices - k) + 1;
            plot(SyllFeat{i}(CorrespondingFirstMotifIndices, j), SyllFeat{i}(Indices,j), 'ko');
            axis tight;
            Temp = axis;
            axis([min(Temp) max(Temp) min(Temp) max(Temp)]);
            plot([min(Temp) max(Temp)], [min(Temp) max(Temp)], 'k--');
            
            MeanSyllFeats{i}(k) = mean(SyllRatioToFirstinBout{i}(Indices, j));
            STDSyllFeats{i}(k) = std(SyllRatioToFirstinBout{i}(Indices, j));
            SEMSyllFeats{i}(k) = std(SyllRatioToFirstinBout{i}(Indices, j))/sqrt(length(FirstMotifSylls));
            
        end
        subplot(PlotNumRows, PlotNumCols, ((j-1)*PlotNumCols) + (k));
        hold on;
        %MeanBar{j} = bar(MeanSyllFeats{i});
        %set(MeanBar{j}, 'FaceColor', 'none', 'EdgeColor', 'k');
        errorbar(MeanSyllFeats{i}, STDSyllFeats{i}, 'ko-');
        plot([0 (k+1)], [1 1], 'k--');
        axis tight;
        Temp = axis;
        axis([0.75 (k+0.25) (Temp(3)*0.98) (Temp(4)*1.02)]);
    end
    
end
    

