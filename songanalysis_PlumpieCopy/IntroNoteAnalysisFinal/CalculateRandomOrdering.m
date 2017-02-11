function [RandomOrderSequences] = CalculateRandomOrdering(AllINFeatLabels, AllINFeats, AllINFeatNoofINs, ValidIntroNoteNos, NumReps)

LastSeqIndex = 1;
FirstSeqIndex = 1;
RandomOrderSequences.CommonLast.SimilarPositionDistances = [];
RandomOrderSequences.CommonLast.OtherPositionDistances = [];
RandomOrderSequences.CommonFirst.SimilarPositionDistances = [];
RandomOrderSequences.CommonFirst.OtherPositionDistances = [];
    
for RepNo = 1:NumReps,
    TempINFeats = [];
    TempINFeatLabels = [];
    TempINFeatNoofINs = [];
    for j = ValidIntroNoteNos',
        Indices = find(AllINFeatNoofINs == j);
        TempIndices = randperm(length(Indices));
        TempINFeatNoofINs = [TempINFeatNoofINs; AllINFeatNoofINs(Indices(TempIndices))];
        TempINFeats = [TempINFeats; AllINFeats(Indices(TempIndices),:)];
        for k = 1:j,
            INGroupIndices = Indices(TempIndices(((k-1)*length(TempIndices)/j + 1):(k*length(TempIndices)/j)));
            TempINFeatLabels = [TempINFeatLabels; [ones(length(INGroupIndices),1)*(k-j) ones(length(INGroupIndices),1)*k]];
        end
    end
    
    INNo = 0;
    for i = ValidIntroNoteNos',
        Indices1 = find(TempINFeatNoofINs == i);
        for j = min(TempINFeatLabels(Indices1,1)):1:max(TempINFeatLabels(Indices1,1)),
            TempSamePosDist = [];
            TempOtherPosDist = [];
            INGroup1 = find((TempINFeatLabels(:,1) == j) & (TempINFeatNoofINs == i));
            for k = ValidIntroNoteNos',
                Indices2 = find(TempINFeatNoofINs == k);
                for m = min(TempINFeatLabels(Indices2,1)):1:max(TempINFeatLabels(Indices2,1)),
                    INGroup2 = find((TempINFeatLabels(:,1) == m) & (TempINFeatNoofINs == k));
                    if ((j == m) && (i == k))
                        continue;
                    else
                        Distance = abs(median(TempINFeats(INGroup1)) - median(TempINFeats(INGroup2)));
                        if (j == m)
                            TempSamePosDist = [TempSamePosDist; Distance];
                            RandomOrderSequences.CommonLast.SimilarPositionDistances = [RandomOrderSequences.CommonLast.SimilarPositionDistances; Distance];
                        else
                            TempOtherPosDist = [TempOtherPosDist; Distance];
                            RandomOrderSequences.CommonLast.OtherPositionDistances = [RandomOrderSequences.CommonLast.OtherPositionDistances; Distance];
                        end
                    end
                end
            end
            if (~isempty(TempSamePosDist) && ~isempty(TempOtherPosDist))
                INNo = INNo + 1;
                if ((length(TempSamePosDist) > 1) && (length(TempOtherPosDist) > 1))
                    RandomOrderSequences.CommonLastIndividualDistances{INNo}(RepNo, :) = [mean(TempSamePosDist) mean(TempOtherPosDist) std(TempSamePosDist) std(TempOtherPosDist) i j];
                else
                    if (length(TempSamePosDist) > 1)
                        RandomOrderSequences.CommonLastIndividualDistances{INNo}(RepNo, :) = [mean(TempSamePosDist) mean(TempOtherPosDist) std(TempSamePosDist) NaN i j];
                    else
                        if (length(TempOtherPosDist) > 1)
                            RandomOrderSequences.CommonLastIndividualDistances{INNo}(RepNo, :) = [mean(TempSamePosDist) mean(TempOtherPosDist) NaN std(TempOtherPosDist) i j];
                        else
                            RandomOrderSequences.CommonLastIndividualDistances{INNo}(RepNo, :) = [mean(TempSamePosDist) mean(TempOtherPosDist) NaN NaN i j];
                        end
                    end
                end
            end
        end
    end
    
    INNo = 0;
    for i = ValidIntroNoteNos',
        Indices1 = find(TempINFeatNoofINs == i);
        for j = min(TempINFeatLabels(Indices1,2)):1:max(TempINFeatLabels(Indices1,2)),
            TempSamePosDist = [];
            TempOtherPosDist = [];
            INGroup1 = find((TempINFeatLabels(:,2) == j) & (TempINFeatNoofINs == i));
            for k = ValidIntroNoteNos',
                Indices2 = find(TempINFeatNoofINs == k);
                for m = min(TempINFeatLabels(Indices2,2)):1:max(TempINFeatLabels(Indices2,2)),
                    INGroup2 = find((TempINFeatLabels(:,2) == m) & (TempINFeatNoofINs == k));
                    if ((j == m) && (i == k))
                        continue;
                    else
                        Distance = abs(median(TempINFeats(INGroup1)) - median(TempINFeats(INGroup2)));
                        if (j == m)
                            TempSamePosDist = [TempSamePosDist; Distance];
                            RandomOrderSequences.CommonFirst.SimilarPositionDistances = [RandomOrderSequences.CommonFirst.SimilarPositionDistances; Distance];
                        else
                            TempOtherPosDist = [TempOtherPosDist; Distance];
                            RandomOrderSequences.CommonFirst.OtherPositionDistances = [RandomOrderSequences.CommonFirst.OtherPositionDistances; Distance];
                        end
                    end
                end
            end
            if (~isempty(TempSamePosDist) && ~isempty(TempOtherPosDist))
                INNo = INNo + 1;
                if ((length(TempSamePosDist) > 1) && (length(TempOtherPosDist) > 1))
                    RandomOrderSequences.CommonFirstIndividualDistances{INNo}(RepNo, :) = [mean(TempSamePosDist) mean(TempOtherPosDist) std(TempSamePosDist) std(TempOtherPosDist) i j];
                else
                    if (length(TempSamePosDist) > 1)
                        RandomOrderSequences.CommonFirstIndividualDistances{INNo}(RepNo, :) = [mean(TempSamePosDist) mean(TempOtherPosDist) std(TempSamePosDist) NaN i j];
                    else
                        if (length(TempOtherPosDist) > 1)
                            RandomOrderSequences.CommonFirstIndividualDistances{INNo}(RepNo, :) = [mean(TempSamePosDist) mean(TempOtherPosDist) NaN std(TempOtherPosDist) i j];
                        else
                            RandomOrderSequences.CommonFirstIndividualDistances{INNo}(RepNo, :) = [mean(TempSamePosDist) mean(TempOtherPosDist) NaN NaN i j];
                        end
                    end
                end
            end
        end
    end
end
disp('Finished random ordering');
