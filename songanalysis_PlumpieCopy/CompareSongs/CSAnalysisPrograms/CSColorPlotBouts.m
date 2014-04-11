function [] = CSColorPlotBouts(CSData)

for i = 1:length(CSData.AllLabels),
    figure; hold on;

    RowIndex = 1;
    
    BoutStarts = find(CSData.AllLabels{i} == 'Q');
    BoutEnds = find(CSData.AllLabels{i} == 'q');
    
    clear FirstMotifSyll NoofINs LastInterval;
    
    for j = 1:length(BoutStarts),
        Flag = 0;
        for k = BoutStarts(j)+1:BoutEnds(j)-1,
            if (~isempty(find(CSData.MotifSyllLabels == CSData.AllLabels{i}(k))))
                Flag = 1;
                break;
            end
        end

        if (Flag == 1)
            FirstMotifSyll(j) = k;
        else
            FirstMotifSyll(j) = 0;
        end

        if (FirstMotifSyll(j) > 1)
            NoofINs(j) = FirstMotifSyll(j) - (BoutStarts(j) + 1);
            LastInterval(j) = CSData.AllOnsets{i}(k) - CSData.AllOffsets{i}(k-1);
        else
            NoofINs(j) = 0;
            LastInterval(j) = 0;
        end
    end

    [SortedData, SortedIndices] = sort(NoofINs);

    for j = SortedIndices,
        if (FirstMotifSyll(j) > 0)

            Onsets = CSData.AllOnsets{i}(BoutStarts(j)+1:BoutEnds(j)-1)/1000;
            Offsets = CSData.AllOffsets{i}(BoutStarts(j)+1:BoutEnds(j)-1)/1000;
            Labels = CSData.AllLabels{i}(BoutStarts(j)+1:BoutEnds(j)-1);
            
            Offsets = Offsets - CSData.AllOnsets{i}(FirstMotifSyll(j))/1000;
            Onsets = Onsets - CSData.AllOnsets{i}(FirstMotifSyll(j))/1000;

            for k = 1:length(Labels),
                if (~isempty(find(CSData.MotifSyllLabels == Labels(k))))
                    plot([Onsets(k) Offsets(k)], [RowIndex RowIndex], 'Color', CSData.MotifLabelsColor, 'LineWidth', 6);
                else
                    if (~isempty(find(CSData.INLabels == Labels(k))))
                        plot([Onsets(k) Offsets(k)], [RowIndex RowIndex], 'Color', CSData.INLabelsColor, 'LineWidth', 6);
                    else
                        plot([Onsets(k) Offsets(k)], [RowIndex RowIndex], 'Color', CSData.OtherSyllsColor, 'LineWidth', 6);
                    end
                end
                if (k ~= length(Labels))
                    plot([Offsets(k) Onsets(k+1)], [RowIndex RowIndex], 'Color', CSData.SilenceColor, 'LineWidth', 6);
                end
            end
            RowIndex = RowIndex + 1;
        end
    end
end

disp('Finished plotting color plots');