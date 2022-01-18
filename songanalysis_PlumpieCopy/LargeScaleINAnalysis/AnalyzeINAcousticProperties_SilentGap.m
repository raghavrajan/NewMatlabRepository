function [] = AnalyzeINAcousticProperties_SilentGap(BirdDetailsTextFile, MinTrialNo)

[BirdParameters] = ProcessSongData(BirdDetailsTextFile, 2000);

% One more pre-processing step. I currently have a list of bouts and I also
% have a separate list of syllables that is compiled across all files. I
% have to have a set of indices that allow me to know which bout each
% syllable in the list of syllables belongs to
for i = 1:length(BirdParameters),
    BirdParameters(i).SyllableListBoutNum = zeros(size(BirdParameters(i).SyllableData,1), 1);
    if (BirdParameters(i).Continuousdata == 0)
        TotalSyllNo = 0;
        for j = 1:length(BirdParameters(i).SongFileNames),
            BoutIndices = find(BirdParameters(i).Bouts(:,3) == j);
            for k = 1:length(BoutIndices),
                BirdParameters(i).SyllableListBoutNum((TotalSyllNo + BirdParameters(i).Bouts(BoutIndices(k), 1)):(TotalSyllNo + BirdParameters(i).Bouts(BoutIndices(k), 2))) = BoutIndices(k);
            end
            TotalSyllNo = TotalSyllNo + length(BirdParameters(i).NoteInfo{j}.labels);
        end
    else
        for j = 1:size(BirdParameters(i).Bouts,1),
            BirdParameters(i).SyllableListBoutNum(BirdParameters(i).Bouts(j,1):BirdParameters(i).Bouts(j,2)) = j;
        end
    end
end

% Now for all the valid bouts - i.e. bouts with atleast 1 motif syllable,
% I'm going to pick out sequences of 2 more INs. Then look at the change in
% properties of the consecutive INs as a function of the interval between
% the two.

for i = 1:length(BirdParameters),
    CurrentIN_SAPFeats{i} = [];
    NextIN_SAPFeats{i} = [];
    ConsecutiveIN_Gaps{i} = [];
    CurrentIN_Position{i} = [];
    
    ValidBouts = find((BirdParameters(i).Bouts(:,7) == 1) & (BirdParameters(i).Bouts(:,8) > 0) & (BirdParameters(i).Bouts(:,9) > 1));
    for j = ValidBouts(:)',
        BoutSyllableIndices = find(BirdParameters(i).SyllableListBoutNum == j);
        BoutLabels = char(BirdParameters(i).SyllableData(BoutSyllableIndices, 1));
        BoutOnsets = BirdParameters(i).SyllableData(BoutSyllableIndices, 4);
        BoutOffsets = BirdParameters(i).SyllableData(BoutSyllableIndices, 5);
        BoutSAPFeats = BirdParameters(i).SAPFeatsMatrix(BoutSyllableIndices,:);
        
        % I also want to keep track of position of the IN in a bout
        % relative to the first motif syllable in the bout. The way I'm
        % going to do this is to keep track of the first motif syllable and
        % then put the positions of the
        MotifSylls = zeros(size(BoutLabels));
        for k = 1:length(BoutLabels),
            if (~isempty(find(BirdParameters(i).MotifLabels == BoutLabels(k))))
                MotifSylls(k) = 1;
            end
        end
        
        
        INs = ones(size(BoutLabels))*NaN;
        for k = 1:length(BoutLabels),
            if (~isempty(find(BirdParameters(i).INLabels == BoutLabels(k))))
                INs(k) = 1;
            end
        end
        
        MotifOnsets = find(conv(MotifSylls, [1 -1]) == 1);
        MotifOffsets = find(conv(MotifSylls, [1 -1]) == -1) - 1;
        
        INPosition = ones(size(BoutLabels))*NaN;
        
        for k = 1:length(MotifOnsets),
            if (k == 1)
                INPosition(1:MotifOnsets(k)-1) = -(MotifOnsets(k) - 1):1:-1;
            else
                INPosition(MotifOffsets(k-1)+1:MotifOnsets(k) - 1) = -(MotifOnsets(k) - MotifOffsets(k-1) - 1):1:-1;
            end
        end
        % Now find consecutive INs
        ConsecutiveINs = find(diff(INs) == 0);
        for k = 1:length(ConsecutiveINs),
            ConsecutiveIN_Gaps{i}(end+1) = BoutOnsets(ConsecutiveINs(k) + 1) - BoutOffsets(ConsecutiveINs(k));
            CurrentIN_SAPFeats{i}(end+1,:) = BoutSAPFeats(ConsecutiveINs(k),:);
            NextIN_SAPFeats{i}(end+1,:) = BoutSAPFeats(ConsecutiveINs(k) + 1,:);
            CurrentIN_Position{i}(end+1,:) = INPosition(k);
        end
    end
    % Now plot the acoustic features as a function of gap size
    figure;
    for j = 1:4,
        subplot(2, 2, j);
        plot(ConsecutiveIN_Gaps{i}, NextIN_SAPFeats{i}(:,j)./CurrentIN_SAPFeats{i}(:,j), 'ko');
        if (j > (length(BirdParameters(i).SAPFeat_FieldNames) - 2))
            xlabel('Gap between INs (msec)');
        end
        ylabel(BirdParameters(i).SAPFeat_FieldNames{j});
        [r, p] = corrcoef(ConsecutiveIN_Gaps{i}, NextIN_SAPFeats{i}(:,j)./CurrentIN_SAPFeats{i}(:,j), 'rows', 'complete');
        if (j < 3)
            if (p(1,2) < 0.05)
                title([BirdParameters(i).BirdName, ':', BirdParameters(i).DataLabel, ':r=', num2str(r(1,2)), ';p=', num2str(p(1,2))]);
            else
                title([BirdParameters(i).BirdName, ':', BirdParameters(i).DataLabel, ':p > 0.05']);
            end
        else
            if (p(1,2) < 0.05)
                title(['r=', num2str(r(1,2)), ';p=', num2str(p(1,2))]);
            else
                title(['p > 0.05']);
            end
        end
    end
        
end


disp('Finished Analysis');