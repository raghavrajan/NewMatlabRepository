function [PC] = IntroNotePopulationNeuralActivityAnalysis(Neural_INR, BinSize, PC, Motif)

Colours = ['rgbcmk'];
PreTime = 0.15;
PostTime = 0.15;

% First get the lengths of all motifs for each site (SU and MU) and get the
% associated lengths of the syllables and gaps

Index = 0;
for i = 1:length(Neural_INR),
    SiteIndex = 0;
    for j = 1:length(Neural_INR{i}.BoutDetails),
        Matches = strfind(Neural_INR{i}.BoutDetails(j).labels, Motif);
        SpikeTimes = Neural_INR{i}.BoutDetails(j).SpikeTimes;
        for k = 1:length(Matches),
            Index = Index + 1;
            SiteIndex = SiteIndex + 1;
            MotifDur{i}(SiteIndex) = Neural_INR{i}.BoutDetails(j).offsets(Matches(k) - 1 + length(Motif)) - Neural_INR{i}.BoutDetails(j).onsets(Matches(k));
            AllMotifDur(Index) = MotifDur{i}(SiteIndex);
            for Sylls = 1:length(Motif),
                TempSyllDur(1, (Sylls - 1)*2 + 1) = Neural_INR{i}.BoutDetails(j).offsets(Matches(k) + Sylls - 1) - Neural_INR{i}.BoutDetails(j).onsets(Matches(k) + Sylls - 1);
            end
            for Gaps = 1:(length(Motif)-1),
                TempSyllDur(1, (Gaps - 1)*2 + 2) = Neural_INR{i}.BoutDetails(j).onsets(Matches(k) + Gaps) - Neural_INR{i}.BoutDetails(j).offsets(Matches(k) + Gaps - 1);
            end
            SyllGapDurs{i}(SiteIndex,:) = TempSyllDur;
            AllSyllGapDurs(Index,:) = TempSyllDur;
        end
    end
end

[SortedAllMotifDurs, SortedIndices] = sort(AllMotifDur);
SortedAllSyllGapDurs = AllSyllGapDurs(SortedIndices,:);

MedianMotif = round(length(SortedIndices)/2);
MedianMotifDur = SortedAllMotifDurs(MedianMotif);
MedianSyllGapDurs = SortedAllSyllGapDurs(MedianMotif, :);
MedianSyllGapOnsets = [0 cumsum(MedianSyllGapDurs)];

Edges = -PreTime:BinSize:(MedianMotifDur + PostTime);

% Now get the spike times for each of the motifs and then warp it to fit
% the reference motif durations

Index = 0;
for i = 1:length(Neural_INR),
    Raster{i} = [];
    SiteIndex = 0;
    for j = 1:length(Neural_INR{i}.BoutDetails),
        Matches = strfind(Neural_INR{i}.BoutDetails(j).labels, Motif);
        SpikeTimes = Neural_INR{i}.BoutDetails(j).SpikeTimes;
        for k = 1:length(Matches),
            Index = Index + 1;
            SiteIndex = SiteIndex + 1;
            TempWarpedSpikeTrain = [];
            MotifOnset = Neural_INR{i}.BoutDetails(j).onsets(Matches(k));
            
            % First get all the spikes associated with the pre-time - there
            % is no warping to be done here
            TempSpikeTimes = SpikeTimes(find((SpikeTimes >= (MotifOnset - PreTime)) & (SpikeTimes < MotifOnset))) - MotifOnset;
            TempSpikeTimes = TempSpikeTimes(:);
            TempWarpedSpikeTrain = [TempWarpedSpikeTrain; TempSpikeTimes];
            
            % Now get the spikes associated with each syllable and each gap
            % and warp it to the correct length
            
            for Sylls = 1:length(Motif),
                % First get the spikes associated with a syllable
                TempSpikeTimes = SpikeTimes(find((SpikeTimes >= Neural_INR{i}.BoutDetails(j).onsets(Matches(k) - 1 + Sylls)) & (SpikeTimes < Neural_INR{i}.BoutDetails(j).offsets(Matches(k) - 1 + Sylls)))) - Neural_INR{i}.BoutDetails(j).onsets(Matches(k) - 1 + Sylls);
                if (~isempty(TempSpikeTimes))
                    TempSpikeTimes = TempSpikeTimes(:);
                    % Now warp the spikes to match the duration of the
                    % reference syllable
                    TempSpikeTimes = TempSpikeTimes * MedianSyllGapDurs(1, (Sylls - 1)*2 + 1)/SyllGapDurs{i}(SiteIndex, (Sylls - 1)*2 + 1);

                    % Now add the onset time of the syllable from the reference
                    % motif to put it in the correct place in the reference
                    % motif
                    TempSpikeTimes = TempSpikeTimes + MedianSyllGapOnsets(1, (Sylls - 1)*2 + 1);

                    TempWarpedSpikeTrain = [TempWarpedSpikeTrain; TempSpikeTimes];
                end             
                % Now get the spikes associated with the next gap
                if (Sylls ~= length(Motif))
                    TempSpikeTimes = SpikeTimes(find((SpikeTimes >= Neural_INR{i}.BoutDetails(j).offsets(Matches(k) - 1 + Sylls)) & (SpikeTimes < Neural_INR{i}.BoutDetails(j).onsets(Matches(k) + Sylls)))) - Neural_INR{i}.BoutDetails(j).offsets(Matches(k) - 1 + Sylls);
                    if (~isempty(TempSpikeTimes))
                        TempSpikeTimes = TempSpikeTimes(:);
                        % Now warp the spikes to match the duration of the
                        % reference gap
                        TempSpikeTimes = TempSpikeTimes * MedianSyllGapDurs(1, (Sylls - 1)*2 + 2)/SyllGapDurs{i}(SiteIndex, (Sylls - 1)*2 + 2);

                        % Now add the onset time of the gap from the reference
                        % motif to put it in the correct place in the reference
                        % motif
                        TempSpikeTimes = TempSpikeTimes + MedianSyllGapOnsets(1, (Sylls - 1)*2 + 2);

                        TempWarpedSpikeTrain = [TempWarpedSpikeTrain; TempSpikeTimes];
                    end
                end
            end
            
            % Now get all the spikes associated with the post-time - there
            % is no warping to be done here
            MotifOffset = Neural_INR{i}.BoutDetails(j).offsets(Matches(k) - 1 + length(Motif));            
            TempSpikeTimes = SpikeTimes(find((SpikeTimes >= MotifOffset) & (SpikeTimes < (MotifOffset + PostTime)))) - MotifOffset + MedianSyllGapOnsets(end);
            TempSpikeTimes = TempSpikeTimes(:);
            TempWarpedSpikeTrain = [TempWarpedSpikeTrain; TempSpikeTimes];
            TempWarpedSpikeTrain = TempWarpedSpikeTrain(:)';
            
            WarpedSpikeTrain{i}{SiteIndex} = TempWarpedSpikeTrain;
            PST{i}(SiteIndex,:) = histc(TempWarpedSpikeTrain, Edges)/BinSize;
            if (~isempty(TempWarpedSpikeTrain))
                Raster{i} = [Raster{i}; [TempWarpedSpikeTrain' ones(size(TempWarpedSpikeTrain'))*SiteIndex]];
            end
        end
    end
    MeanPST(i,:) = mean(PST{i}(:,1:end-1))/std(PST{i}(1:size(PST{i},1)*size(PST{i},2)));
end
   
if (isempty(PC))
    [PC, Var] = eig(cov(MeanPST'));
end

MeanScore = MeanPST'*PC;

figure;
hold on;
for i = 1:length(MedianSyllGapOnsets),
    Index = find(Edges >= MedianSyllGapOnsets(i), 1, 'first');
    if (i == 1)
        Segments(i,:) = [1 Index];
    else
        Segments(i,:) = [Segments(i-1,2) Index];
    end
end
Segments(i+1,:) = [Segments(i,2) length(Edges)-1];

for j = 1:1,
    clear TempPST;
    for i = 1:length(PST),
        TempPST(i,:) = mean(PST{i}(j:j+2,1:end-1))/std(PST{i}(1:size(PST{i},1)*size(PST{i},2)));
    end
    Score = TempPST'*PC;
    
    SyllIndex = 0;
    for i = 1:size(Segments,1),
        if (i == 1)
            PlotCol = [Colours(1), '-'];
        else
            if (i == size(Segments,1))
                PlotCol = [Colours(1), '-'];
            else
                if (mod(i,2) == 0)
                    SyllIndex = SyllIndex + 1;
                    PlotCol = [Colours(1+SyllIndex), '-'];
                else
                    PlotCol = [Colours(end), '-'];
                end
            end
        end
        plot3(Score(Segments(i,1):Segments(i,2), end), Score(Segments(i,1):Segments(i,2), end-1) , Score(Segments(i,1):Segments(i,2), end-2), 'Color', [0.75 0.75 0.75]);
        if (j == 5)
            plot3(MeanScore(Segments(i,1):Segments(i,2), end), MeanScore(Segments(i,1):Segments(i,2), end-1) , MeanScore(Segments(i,1):Segments(i,2), end-2), PlotCol, 'LineWidth', 2);
        end
    end
end
Score = MeanPST'*PC;
disp('Finished');
