function [] = VishrutaShikha_PlotMeanINNumbers_TutorvsPupil(BirdParameters)


for i = 1:length(BirdParameters),
    % First find valid song bouts
    disp(['Bird #', num2str(i)]);
    SongBouts = find((BirdParameters(i).Bouts(:,7) == 1) & (BirdParameters(i).Bouts(:,8) >= 0) & (BirdParameters(i).Bouts(:,9) >= 1));
    
    % now for each of these bouts, find the number of INs at the beginning.
    fprintf('\n');
    for j = 1:length(SongBouts),
        fprintf('%d > ', j);
        % First find first motif syllable
        BoutIndices = find(BirdParameters(i).SyllableListBoutNum == SongBouts(j));
        BoutLabels = char(BirdParameters(i).SyllableData(BoutIndices,1));
        
        for k = 1:length(BoutLabels),
            if (~isempty(find(BirdParameters(i).MotifLabels == BoutLabels(k))))
                OverallFirstMotifSyllIndex = k + (BoutIndices(1) - 1);
                FirstMotifSyllIndex = k;
                break;
            end
        end
        
        if (BirdParameters(i).Continuousdata == 1)
            MotifOnsetTime{i}(j) = BirdParameters(i).SyllableData(OverallFirstMotifSyllIndex,6) - BirdParameters(i).SyllableData(BoutIndices(1), 6);
        else
            MotifOnsetTime{i}(j) = BirdParameters(i).SyllableData(OverallFirstMotifSyllIndex,4) - BirdParameters(i).SyllableData(BoutIndices(1), 4);
        end
        
        INs{i}{j} = zeros(FirstMotifSyllIndex - 1, 1);
        for k = 1:length(BoutLabels(1:(FirstMotifSyllIndex - 1)))
            if (~isempty(find(BirdParameters(i).INLabels == BoutLabels(k))))
                INs{i}{j}(k) = 1;
            end
        end
    
        % Now we have to calculate # of INs for each bout
        % We'll do it two ways
        % a) just total # of INs in each bout
        % b) the last set of consecutive INs that have < 500ms gap between them
    
        TotalNumINs{i}(j) = length(find(INs{i}{j} == 1));
        
        % First to ensure that we only take the last set of consecutive INs
        if (TotalNumINs{i}(j) > 0)
            Gaps = diff(INs{i}{j});
            LongGaps = find(Gaps > 1);
            if (~isempty(LongGaps))
                ConsecutiveINIndices = cumsum(INs{i}{j}(LongGaps(end)+1:end)) + LongGaps(end);
            else
                ConsecutiveINIndices = cumsum(INs{i}{j});
            end

            % Now to check that the gap between these INs is not more than
            % 500ms

            if (BirdParameters(i).Continuousdata == 1)
                Gaps = BirdParameters(i).SyllableData((ConsecutiveINIndices(2:end) + BoutIndices(1)), 6) - BirdParameters(i).SyllableData((ConsecutiveINIndices(1:end-1) + BoutIndices(1)), 7);
            else 
                Gaps = BirdParameters(i).SyllableData((ConsecutiveINIndices(2:end) + BoutIndices(1)), 4) - BirdParameters(i).SyllableData((ConsecutiveINIndices(1:end-1) + BoutIndices(1)), 5);
            end

            LongGaps = find(Gaps > 500);
            if (~isempty(LongGaps))
                ConsecutiveINIndices = ConsecutiveINIndices(LongGaps(end)+1:end);
            end
            NumConsecutiveINs{i}(j) = length(ConsecutiveINIndices);
        end    
    end
end



disp('Finished');

