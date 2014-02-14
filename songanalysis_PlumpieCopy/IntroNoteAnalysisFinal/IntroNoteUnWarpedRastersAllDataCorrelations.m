function [PairWiseCorr, RandomPairWiseCorr, LastINCorr, RandomLastINCorr, RandomLastINCorr2, Position, StartTimes, Edges, INDur, GapDur] = IntroNoteUnWarpedRastersAllDataCorrelations(Neural_INR, PreTime, PostTime, MinPlotNumber, varargin)

if (nargin > 4)
    if (~isempty(varargin{1}))
        PlotOption = varargin{1};
    else
        PlotOption = 'on';
    end
    if (nargin > 5)
        ComparisonIN = varargin{2};
    else
        ComparisonIN = -1;
    end
else
    PlotOption = 'on';
    ComparisonIN = -1;
end

Fs = 10000;
Width = 0.01;
GaussianLen = 4;

XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

Edges = -PreTime:1/Fs:PostTime;

INIndex = 0;
TotalTrials = sum(Neural_INR.NoofINs(find(Neural_INR.NoofINs > 0))) + sum(Neural_INR.WithinBoutNoofINs(find(Neural_INR.WithinBoutNoofINs(:,1) > 0)));
SpikeTrain = zeros(TotalTrials, length(Edges));

for i = 1:length(Neural_INR.NoofINs),
    if (Neural_INR.NoofINs(i) > 0)
        INs = Neural_INR.INs{i};
        for j = 1:length(INs),
            INIndex = INIndex + 1;
            
            Position(INIndex,1) = j - length(INs) - 1;
            Position(INIndex,2) = j;
            NoofINs(INIndex,1) = length(INs);
            
            INOnset = Neural_INR.BoutDetails(i).onsets(INs(j));
            INOffset = Neural_INR.BoutDetails(i).offsets(INs(j));
            GapOnset = INOffset;
            GapOffset = Neural_INR.BoutDetails(i).onsets(INs(j)+1);
            
            INDur(INIndex) = INOffset - INOnset;
            INGapDur(INIndex) = GapOffset - INOnset;
            
            %StartIndex = find(Time{i} >= (INOnset + Edges(1)), 1, 'first');
            
            SpikeTimes{INIndex} = Neural_INR.BoutDetails(i).SpikeTimes(find((Neural_INR.BoutDetails(i).SpikeTimes >= (INOnset + Edges(1))) & (Neural_INR.BoutDetails(i).SpikeTimes < (INOnset + Edges(end))))) - INOnset;
            SpikeTrainIndices = round((SpikeTimes{INIndex} - Edges(1)) * Fs);
            SpikeTrainIndices(find(SpikeTrainIndices == 0)) = 1;
            SpikeTrain(INIndex, SpikeTrainIndices) = 1;
            SpikeTrain(INIndex,:) = conv(SpikeTrain(INIndex,:), GaussWin, 'same');
        end
    end
end

for i = 1:size(Neural_INR.WithinBoutNoofINs, 1),
    if (Neural_INR.WithinBoutNoofINs(i,1) > 0)
        BoutSpikeTimes = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes;
        INs = Neural_INR.WithinBoutINs{i};
        for j = 1:length(INs),
            INIndex = INIndex + 1;
            Position(INIndex,1) = j - length(INs) - 1;
            Position(INIndex,2) = j;
            NoofINs(INIndex,1) = length(INs);
            
            INOnset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j));
            INOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(INs(j));
            GapOnset = INOffset;
            GapOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j)+1);
            
            INDur(INIndex) = INOffset - INOnset;
            INGapDur(INIndex) = GapOffset - INOnset;
            
            SpikeTimes{INIndex} = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes(find((Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes >= (INOnset + Edges(1))) & (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes < (INOnset + Edges(end))))) - INOnset;
            SpikeTrainIndices = round((SpikeTimes{INIndex} - Edges(1)) * Fs);
            SpikeTrainIndices(find(SpikeTrainIndices == 0)) = 1;
            SpikeTrain(INIndex, SpikeTrainIndices) = 1;
            SpikeTrain(INIndex,:) = conv(SpikeTrain(INIndex,:), GaussWin, 'same');
            
            %StartIndex = find(Time{Neural_INR.WithinBoutINBoutIndices(i)} >= (INOnset + Edges(1)), 1, 'first');
            
            %SpikeTrain(INIndex, :) = FR{Neural_INR.WithinBoutINBoutIndices(i)}(StartIndex:(StartIndex - 1 + length(Edges)));

        end
    end
end

% First calculate pairwise correlations between all pairs of spike trains
% for INs at a given position

for i = min(Position(:,1)):1:max(Position(:,1)),
    if (length(find(Position(:,1) == i)) < 3)
        continue;
    end
    SpecificINSpikeTrains = SpikeTrain(find(Position(:,1) == i),:) - repmat(mean(SpikeTrain(find(Position(:,1) == i),:),2), 1, size(SpikeTrain,2));
    clear NormINSpikeTrains;
    for j = 1:size(SpecificINSpikeTrains, 1),
        NormINSpikeTrains(j,1) = norm(SpecificINSpikeTrains(j,:));
    end
    NormINSpikeTrains = repmat(NormINSpikeTrains, 1, size(SpecificINSpikeTrains,1));
    
    PairWiseCorr{abs(i)} = SpecificINSpikeTrains*SpecificINSpikeTrains'./((NormINSpikeTrains).*(NormINSpikeTrains'));
    PairWiseCorr{abs(i)} = triu(PairWiseCorr{abs(i)},1);
    PairWiseCorr{abs(i)} = PairWiseCorr{abs(i)}(:);
    PairWiseCorr{abs(i)}(find(PairWiseCorr{abs(i)} == 0)) = [];
    PairWiseCorr{abs(i)}(isnan(PairWiseCorr{abs(i)})) = 0;
end

% Now calculate pairwise correlations between all pairs of spike trains
% for INs at a given position after each spike train has been shifted by a
% random amount

RandomSpikeTrain = zeros(size(SpikeTrain));
for j = 1:length(SpikeTimes),
    TempSpikeTimes{j} = SpikeTimes{j} + (rand * (Edges(end) - Edges(1)));
    TempSpikeTimes{j}(find(TempSpikeTimes{j} > Edges(end))) = TempSpikeTimes{j}(find(TempSpikeTimes{j} > Edges(end))) - (Edges(end) - Edges(1));
    SpikeTrainIndices = round((TempSpikeTimes{j} - Edges(1)) * Fs);
    SpikeTrainIndices(find(SpikeTrainIndices == 0)) = 1;
    RandomSpikeTrain(j, SpikeTrainIndices) = 1;
    RandomSpikeTrain(j,:) = conv(RandomSpikeTrain(j,:), GaussWin, 'same');
end

for i = min(Position(:,1)):1:max(Position(:,1)),
    if (length(find(Position(:,1) == i)) < 3)
        continue;
    end
      
    SpecificINSpikeTrains = RandomSpikeTrain(find(Position(:,1) == i),:) - repmat(mean(RandomSpikeTrain(find(Position(:,1) == i),:),2), 1, size(RandomSpikeTrain,2));
    clear NormINSpikeTrains;
    for j = 1:size(SpecificINSpikeTrains, 1),
        NormINSpikeTrains(j,1) = norm(SpecificINSpikeTrains(j,:));
    end
    NormINSpikeTrains = repmat(NormINSpikeTrains, 1, size(SpecificINSpikeTrains,1));
    
    RandomPairWiseCorr{abs(i)} = SpecificINSpikeTrains*SpecificINSpikeTrains'./((NormINSpikeTrains).*(NormINSpikeTrains'));
    RandomPairWiseCorr{abs(i)} = triu(RandomPairWiseCorr{abs(i)},1);
    RandomPairWiseCorr{abs(i)} = RandomPairWiseCorr{abs(i)}(:);
    RandomPairWiseCorr{abs(i)}(find(RandomPairWiseCorr{abs(i)} == 0)) = [];
    RandomPairWiseCorr{abs(i)}(isnan(RandomPairWiseCorr{abs(i)})) = 0;
end

BinSize = 0.04;

Edges = -PreTime:1/Fs:(PostTime + BinSize);

INIndex = 0;
TotalTrials = sum(Neural_INR.NoofINs(find(Neural_INR.NoofINs > 0))) + sum(Neural_INR.WithinBoutNoofINs(find(Neural_INR.WithinBoutNoofINs(:,1) > 0)));
SpikeTrain = zeros(TotalTrials, length(Edges));

SyllDur = [];
for i = 1:length(Neural_INR.NoofINs),
    if (Neural_INR.NoofINs(i) > 0)
        INs = Neural_INR.INs{i};
        for j = 1:length(INs),
            INIndex = INIndex + 1;
            
            Position(INIndex,1) = j - length(INs) - 1;
            Position(INIndex,2) = j;
            NoofINs(INIndex,1) = length(INs);
            
            INOnset = Neural_INR.BoutDetails(i).onsets(INs(j));
            INOffset = Neural_INR.BoutDetails(i).offsets(INs(j));
            GapOnset = INOffset;
            GapOffset = Neural_INR.BoutDetails(i).onsets(INs(j)+1);
            
            INDur(INIndex) = INOffset - INOnset;
            GapDur(INIndex) = GapOffset - GapOnset;
            INGapDur(INIndex) = GapOffset - INOnset;

            if (j == length(INs))
                SyllOnset = Neural_INR.BoutDetails(i).onsets(INs(j) + 1);
                SyllOffset = Neural_INR.BoutDetails(i).offsets(INs(j) + 1);
                SyllDur = [SyllDur; (SyllOffset - SyllOnset)];
            end
            %StartIndex = find(Time{i} >= (INOnset + Edges(1)), 1, 'first');
            
            SpikeTimes{INIndex} = Neural_INR.BoutDetails(i).SpikeTimes(find((Neural_INR.BoutDetails(i).SpikeTimes >= (INOnset + Edges(1))) & (Neural_INR.BoutDetails(i).SpikeTimes < (INOnset + Edges(end))))) - INOnset;
            SpikeTrainIndices = round((SpikeTimes{INIndex} - Edges(1)) * Fs);
            SpikeTrainIndices(find(SpikeTrainIndices == 0)) = 1;
            SpikeTrain(INIndex, SpikeTrainIndices) = 1;
            SpikeTrain(INIndex,:) = conv(SpikeTrain(INIndex,:), GaussWin, 'same');
        end
    end
end

for i = 1:size(Neural_INR.WithinBoutNoofINs, 1),
    if (Neural_INR.WithinBoutNoofINs(i,1) > 0)
        BoutSpikeTimes = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes;
        INs = Neural_INR.WithinBoutINs{i};
        for j = 1:length(INs),
            INIndex = INIndex + 1;
            Position(INIndex,1) = j - length(INs) - 1;
            Position(INIndex,2) = j;
            NoofINs(INIndex,1) = length(INs);
            
            INOnset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j));
            INOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(INs(j));
            GapOnset = INOffset;
            GapOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j)+1);
            
            INDur(INIndex) = INOffset - INOnset;
            GapDur(INIndex) = GapOffset - GapOnset;
            INGapDur(INIndex) = GapOffset - INOnset;
            
            if (j == length(INs))
                SyllOnset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j) + 1);
                SyllOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(INs(j) + 1);
                SyllDur = [SyllDur; (SyllOffset - SyllOnset)];
            end
            
            SpikeTimes{INIndex} = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes(find((Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes >= (INOnset + Edges(1))) & (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes < (INOnset + Edges(end))))) - INOnset;
            SpikeTrainIndices = round((SpikeTimes{INIndex} - Edges(1)) * Fs);
            SpikeTrainIndices(find(SpikeTrainIndices == 0)) = 1;
            SpikeTrain(INIndex, SpikeTrainIndices) = 1;
            SpikeTrain(INIndex,:) = conv(SpikeTrain(INIndex,:), GaussWin, 'same');
            
            %StartIndex = find(Time{Neural_INR.WithinBoutINBoutIndices(i)} >= (INOnset + Edges(1)), 1, 'first');
            
            %SpikeTrain(INIndex, :) = FR{Neural_INR.WithinBoutINBoutIndices(i)}(StartIndex:(StartIndex - 1 + length(Edges)));

        end
    end
end

% Two ways of assessing significance
% 1) Choose all the last INs and repeat it such that I have 1000 last INs
% and then add a random time-shift to each one of them and calculate
% correlation - this would represent the correlations between random spike
% trains and mean last IN firing pattern
% 2) If there are 150 last INs, then take 150 INs at random disregarding
% their proper position and calculate the mean correlation to last IN
% firing pattern - do this 1000 times and choose the top 5%ile as the
% cut-off - this would represent the correlations if position didn't matter

StartTimes = 1:round(0.005*Fs):(size(SpikeTrain,2) - round(BinSize*Fs));

for i = 1:size(Position,1),
    clear Indices;
    if (ComparisonIN > 0)
        Indices = (Position(:,2) == ComparisonIN);
        if (ComparisonIN == Position(i,2))
            Indices(i) = 0;
        end
    else
        Indices = (Position(:,1) == ComparisonIN);
        if (ComparisonIN == Position(i,1))
            Indices(i) = 0;
        end
    end
    Indices = find(Indices > 0);
    
    TrialMeanComparisonIN(i,:) = mean(SpikeTrain(Indices,:));
    
    ColIndex = 0;
    for j = 1:length(StartTimes)
        ColIndex = ColIndex + 1;
        MeanComparisonIN{ColIndex}(i,:) = TrialMeanComparisonIN(i, StartTimes(j):(StartTimes(j) + round(BinSize*Fs)));
        MeanComparisonIN{ColIndex}(i,:) = MeanComparisonIN{ColIndex}(i,:) - mean(MeanComparisonIN{ColIndex}(i,:));
        
        NormMeanComparisonIN{ColIndex}(i,1) = norm(MeanComparisonIN{ColIndex}(i,:));
        
        TrialSpikeTrain{ColIndex}(i,:) = SpikeTrain(i, StartTimes(j):(StartTimes(j) + round(BinSize*Fs)));
        TrialSpikeTrain{ColIndex}(i,:) = TrialSpikeTrain{ColIndex}(i,:) - mean(TrialSpikeTrain{ColIndex}(i,:));
        
        NormTrialSpikeTrain{ColIndex}(i,1) = norm(TrialSpikeTrain{ColIndex}(i,:));
    end
end

for i = 1:length(StartTimes),
    LastINCorr(:,i) = diag(TrialSpikeTrain{i}*MeanComparisonIN{i}')./(NormMeanComparisonIN{i}.*NormTrialSpikeTrain{i});
end

ColIndex = 0;
if (ComparisonIN > 0)
    TrialMeanAllComparisonIN(1,:) = mean(SpikeTrain(find(Position(:,2) == ComparisonIN),:));
else
    TrialMeanAllComparisonIN(1,:) = mean(SpikeTrain(find(Position(:,1) == ComparisonIN),:));
end

for j = 1:length(StartTimes)
    ColIndex = ColIndex + 1;
    MeanAllComparisonIN{ColIndex}(1,:) = TrialMeanAllComparisonIN(1, StartTimes(j):(StartTimes(j) + round(BinSize*Fs)));
    MeanAllComparisonIN{ColIndex}(1,:) = MeanAllComparisonIN{ColIndex}(1,:) - mean(MeanAllComparisonIN{ColIndex}(1,:));

    if (ComparisonIN > 0)
        MeanAllComparisonIN{ColIndex} = repmat(MeanAllComparisonIN{ColIndex}, length(find(Position(:,2) == ComparisonIN)), 1);
    else
        MeanAllComparisonIN{ColIndex} = repmat(MeanAllComparisonIN{ColIndex}, length(find(Position(:,1) == ComparisonIN)), 1);
    end

    NormMeanAllComparisonIN{ColIndex}(1,1) = norm(MeanAllComparisonIN{ColIndex}(1,:));
    if (ComparisonIN > 0)
        NormMeanAllComparisonIN{ColIndex} = repmat(NormMeanAllComparisonIN{ColIndex}, length(find(Position(:,2) == ComparisonIN)), 1);
    else
        NormMeanAllComparisonIN{ColIndex} = repmat(NormMeanAllComparisonIN{ColIndex}, length(find(Position(:,1) == ComparisonIN)), 1);
    end
end
    
%RandomLastINCorr = [];
for RandomTrialNo = 1:1000,
    if (ComparisonIN > 0)
        LastINIndices = find(Position(:,2) == ComparisonIN);
    else
        LastINIndices = find(Position(:,1) == ComparisonIN);
    end
    RandomSpikeTrain = zeros(length(LastINIndices), size(SpikeTrain,2));
    RandomOrder = LastINIndices;

    for j = 1:length(RandomOrder),
        TempSpikeTimes{j} = SpikeTimes{RandomOrder(j)} + (rand * (Edges(end) - Edges(1)));
        TempSpikeTimes{j}(find(TempSpikeTimes{j} > Edges(end))) = TempSpikeTimes{j}(find(TempSpikeTimes{j} > Edges(end))) - (Edges(end) - Edges(1));
        SpikeTrainIndices = round((TempSpikeTimes{j} - Edges(1)) * Fs);
        SpikeTrainIndices(find(SpikeTrainIndices == 0)) = 1;
        RandomSpikeTrain(j, SpikeTrainIndices) = 1;
        RandomSpikeTrain(j,:) = conv(RandomSpikeTrain(j,:), GaussWin, 'same');
    end
    
    ColIndex = 0;
    for j = 1:length(StartTimes)
        ColIndex = ColIndex + 1;
        TempTrialSpikeTrain{ColIndex} = RandomSpikeTrain(:, StartTimes(j):(StartTimes(j) + round(BinSize*Fs))); 
        TempTrialSpikeTrain{ColIndex} = TempTrialSpikeTrain{ColIndex} - repmat(mean(TempTrialSpikeTrain{ColIndex}, 2), 1, size(TempTrialSpikeTrain{ColIndex}, 2));
        for k = 1:size(TempTrialSpikeTrain{ColIndex},1),
            NormTempTrialSpikeTrain{ColIndex}(k,1) = norm(TempTrialSpikeTrain{ColIndex}(k,:));
        end
        TempRandomLastINCorr(:,j) = diag(TempTrialSpikeTrain{j}*MeanAllComparisonIN{j}')./(NormMeanAllComparisonIN{j}.*NormTempTrialSpikeTrain{j});
    end
    RandomLastINCorr(RandomTrialNo,:) = mean(TempRandomLastINCorr);
end
RandomLastINCorr(isnan(RandomLastINCorr)) = 0;
RandomLastINCorr = sort(RandomLastINCorr);

RandomLastINCorr2 = [];
% for RandomTrialNo = 1:1000,
%     RandomIndices = randperm(size(Position,1));
%     LastINIndices = find(Position(:,1) == ComparisonIN);
%     RandomSpikeTrain = SpikeTrain(RandomIndices(1:length(LastINIndices)),:);
% 
%     ColIndex = 0;
%     for j = 1:length(StartTimes)
%         ColIndex = ColIndex + 1;
%         TempTrialSpikeTrain{ColIndex} = RandomSpikeTrain(:, StartTimes(j):(StartTimes(j) + round(BinSize*Fs))); 
%         TempTrialSpikeTrain{ColIndex} = TempTrialSpikeTrain{ColIndex} - repmat(mean(TempTrialSpikeTrain{ColIndex}, 2), 1, size(TempTrialSpikeTrain{ColIndex}, 2));
%         for k = 1:size(TempTrialSpikeTrain{ColIndex},1),
%             NormTempTrialSpikeTrain{ColIndex}(k,1) = norm(TempTrialSpikeTrain{ColIndex}(k,:));
%         end
%         TempRandomLastINCorr(:,j) = diag(TempTrialSpikeTrain{j}*MeanAllComparisonIN{j}')./(NormMeanAllComparisonIN{j}.*NormTempTrialSpikeTrain{j});
%     end
%     RandomLastINCorr2(RandomTrialNo,:) = mean(TempRandomLastINCorr);
% end
% RandomLastINCorr2(isnan(RandomLastINCorr2)) = 0;
% RandomLastINCorr2 = sort(RandomLastINCorr2);


for i = min(Position(:,1)):1:max(Position(:,1)),
    TrialNos(abs(i)) = length(find(Position(:,1) == i));
end

LastINCorr(isnan(LastINCorr)) = 0;

if (strfind(PlotOption, 'on'))
    
    CorrFig = figure;
    set(gcf, 'Color', 'w');

    MaxPlots = length(find(TrialNos >= MinPlotNumber));
    
    PlotWGap = 0.035;

    PlotHt = 0.8;

    PlotWidth = 0.96/MaxPlots - PlotWGap;

    Index = 0;
    
    for i = min(Position(:,1)):1:max(Position(:,1)),
        if (TrialNos(abs(i)) >= MinPlotNumber)
            Index = Index + 1;
            subplot('Position', [(1 - (abs(i))*(PlotWidth + PlotWGap)) 0.15 PlotWidth PlotHt]);
            hold on;
            MeanCorr = mean(LastINCorr(find(Position(:,1) == i),:));
            STDCorr = std(LastINCorr(find(Position(:,1) == i),:));
            MeanINDur = median(INDur(find(Position(:,1) == i)));
            fill([0 MeanINDur MeanINDur 0], [1.5 1.5 -1.5 -1.5], 'k', 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
            if (i > min(Position(:,1)))
                PrevINStart = 0 - median(GapDur(find(Position(:,1) == (i-1)))) - median(INDur(find(Position(:,1) == (i-1))));
                PrevINEnd = 0 - median(GapDur(find(Position(:,1) == (i-1))));
                fill([PrevINStart PrevINEnd PrevINEnd PrevINStart], [1.5 1.5 -1.5 -1.5], 'k', 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);                
            end
            
            if (i < max(Position(:,1)))
                NextINStart = MeanINDur + median(GapDur(find(Position(:,1) == i)));
                NextINEnd = MeanINDur + median(GapDur(find(Position(:,1) == i))) + median(INDur(find(Position(:,1) == (i+1))));
                fill([NextINStart NextINEnd NextINEnd NextINStart], [1.5 1.5 -1.5 -1.5], 'k', 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
            else
                NextINStart = MeanINDur + median(GapDur(find(Position(:,1) == i)));
                NextINEnd = MeanINDur + median(GapDur(find(Position(:,1) == i))) + median(SyllDur);
                fill([NextINStart NextINEnd NextINEnd NextINStart], [1.5 1.5 -1.5 -1.5], 'k', 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
            end
            
            fill([(Edges(StartTimes)) fliplr(Edges(StartTimes) + 0.005/2)], [(MeanCorr + STDCorr) fliplr(MeanCorr - STDCorr)], 'k', 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
            plot(Edges(StartTimes), MeanCorr, 'k');
            plot(Edges(StartTimes), RandomLastINCorr(round(0.95*1000),:), 'k--');
            %plot(Edges(StartTimes) + 0.005/2, RandomLastINCorr2(round(0.95*1000),:), 'b--');
            
%            plot([-PreTime PostTime], [1 1], 'k--');
%            plot([-PreTime PostTime], [0 0], 'k--');
%            plot([-PreTime PostTime], [-1 -1], 'k--');
%            set(gca, 'YColor', 'w');
            set(gca, 'Box', 'off');
            axis([-PreTime PostTime -1.255 1.255]);
            set(gca, 'TickDir', 'out');
            set(gca, 'YTick', [-1 0 1], 'YTickLabel', []);
            set(gca, 'XTick', [-PreTime:0.1:PostTime], 'XTickLabel', []);
            set(gca, 'TickLen', [0.05 0.025]);
        end
    end
end
disp('Finished');