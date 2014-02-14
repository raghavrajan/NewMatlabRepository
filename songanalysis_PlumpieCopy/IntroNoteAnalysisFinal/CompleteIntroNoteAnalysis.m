function [Position, Results, PST, SpikeTimes, Stats] = CompleteIntroNoteAnalysis(Neural_INR, BinSize, PreMotorLatency, ComparisonIN)

%=========================CompleteIntroNoteAnalysis Help==================%
% Usage: 
% [Position, Results, Stats] = CompleteIntroNoteAnalysis(Neural_INR, BinSize, PreMotorLatency);

% Position is organised as follows:
% Column 1 - position of intro note relative to last IN
% Column 2 - position of intro note relative to first IN
% Column 3 - no of intro notes in sequence
% Column 4 - 1 for bout beginnings and 0 for within bouts

% Results is organised as follows:
% Column 1 - firing rate during the intro note
% Column 2 - duration of the intro note
% Column 3 - firing rate during the gap
% Column 4 - duration of the gap
% Column 5 - firing rate during the combined intro note and gap
% Column 6 - duration of the combined intro note and gap
% Column 7 - PST during the intro note
% Column 8 - PST during the gap
% Column 9 - PST during the combined intro note and gap
%=========================================================================%

FFWindow = 0.03;
FFAdvance = 0.001;

Width = 0.005;
GaussianLen = 2;
IFRFs = 1/BinSize;
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (IFRFs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * IFRFs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * IFRFs) * (Width * IFRFs)));

INIndex = 0;
for i = 1:length(Neural_INR.NoofINs),
    if (Neural_INR.NoofINs(i) > 0)
        BoutSpikeTimes = Neural_INR.BoutDetails(i).SpikeTimes;
        INs = Neural_INR.INs{i};
        for j = 1:length(INs),
            INIndex = INIndex + 1;
            Position(INIndex, :) = [(j - length(INs) - 1) j length(INs) 1];
            INOnset = Neural_INR.BoutDetails(i).onsets(INs(j)) - PreMotorLatency;
            INOffset = Neural_INR.BoutDetails(i).offsets(INs(j)) - PreMotorLatency;
            GapOnset = INOffset;
            GapOffset = Neural_INR.BoutDetails(i).onsets(INs(j)+1) - PreMotorLatency;
            
            Results(INIndex, 1) = length(find((BoutSpikeTimes >= INOnset) & (BoutSpikeTimes < INOffset)))/(INOffset - INOnset);
            Results(INIndex, 2) = INOffset - INOnset;
            
            Results(INIndex, 3) = length(find((BoutSpikeTimes >= GapOnset) & (BoutSpikeTimes < GapOffset)))/(GapOffset - GapOnset);
            Results(INIndex, 4) = GapOffset - GapOnset;
            
            Results(INIndex, 5) = length(find((BoutSpikeTimes >= INOnset) & (BoutSpikeTimes < GapOffset)))/(GapOffset - INOnset);
            Results(INIndex, 6) = GapOffset - INOnset;

            PST{1}{INIndex} = histc(BoutSpikeTimes, (INOnset:BinSize:INOffset))/BinSize;
            PST{1}{INIndex}(end) = [];
            PST{1}{INIndex} = conv(PST{1}{INIndex}, GaussWin, 'same');
            
            PST{2}{INIndex} = histc(BoutSpikeTimes, (GapOnset:BinSize:GapOffset))/BinSize;
            PST{2}{INIndex}(end) = [];
            PST{2}{INIndex} = conv(PST{2}{INIndex}, GaussWin, 'same');
            
            PST{3}{INIndex} = histc(BoutSpikeTimes, (INOnset:BinSize:GapOffset))/BinSize;
            PST{3}{INIndex}(end) = [];
            PST{3}{INIndex} = conv(PST{3}{INIndex}, GaussWin, 'same');
            
            SpikeTimes{1}{INIndex} = BoutSpikeTimes(find((BoutSpikeTimes >= INOnset) & (BoutSpikeTimes < INOffset))) - INOnset;
            SpikeTimes{2}{INIndex} = BoutSpikeTimes(find((BoutSpikeTimes >= GapOnset) & (BoutSpikeTimes < GapOffset))) - GapOnset;
            SpikeTimes{3}{INIndex} = BoutSpikeTimes(find((BoutSpikeTimes >= INOnset) & (BoutSpikeTimes < GapOffset))) - INOnset;
        end
    end
end

for i = 1:size(Neural_INR.WithinBoutNoofINs, 1),
    if (Neural_INR.WithinBoutNoofINs(i,1) > 0)
        BoutSpikeTimes = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes;
        INs = Neural_INR.WithinBoutINs{i};
        for j = 1:length(INs),
            INIndex = INIndex + 1;
            Position(INIndex, :) = [(j - length(INs) - 1) j length(INs) 0];
            INOnset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j)) - PreMotorLatency;
            INOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(INs(j)) - PreMotorLatency;
            GapOnset = INOffset;
            GapOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j)+1) - PreMotorLatency;
            
            Results(INIndex, 1) = length(find((BoutSpikeTimes >= INOnset) & (BoutSpikeTimes < INOffset)))/(INOffset - INOnset);
            Results(INIndex, 2) = INOffset - INOnset;
            
            Results(INIndex, 3) = length(find((BoutSpikeTimes >= GapOnset) & (BoutSpikeTimes < GapOffset)))/(GapOffset - GapOnset);
            Results(INIndex, 4) = GapOffset - GapOnset;
            
            Results(INIndex, 5) = length(find((BoutSpikeTimes >= INOnset) & (BoutSpikeTimes < GapOffset)))/(GapOffset - INOnset);
            Results(INIndex, 6) = GapOffset - INOnset;

            PST{1}{INIndex} = histc(BoutSpikeTimes, (INOnset:BinSize:INOffset))/BinSize;
            PST{1}{INIndex}(end) = [];
            PST{1}{INIndex} = conv(PST{1}{INIndex}, GaussWin, 'same');
            
            PST{2}{INIndex} = histc(BoutSpikeTimes, (GapOnset:BinSize:GapOffset))/BinSize;
            PST{2}{INIndex}(end) = [];
            PST{2}{INIndex} = conv(PST{2}{INIndex}, GaussWin, 'same');
            
            PST{3}{INIndex} = histc(BoutSpikeTimes, (INOnset:BinSize:GapOffset))/BinSize;
            PST{3}{INIndex}(end) = [];
            PST{3}{INIndex} = conv(PST{3}{INIndex}, GaussWin, 'same');
            
            SpikeTimes{1}{INIndex} = BoutSpikeTimes(find((BoutSpikeTimes >= INOnset) & (BoutSpikeTimes < INOffset))) - INOnset;
            SpikeTimes{2}{INIndex} = BoutSpikeTimes(find((BoutSpikeTimes >= GapOnset) & (BoutSpikeTimes < GapOffset))) - GapOnset;
            SpikeTimes{3}{INIndex} = BoutSpikeTimes(find((BoutSpikeTimes >= INOnset) & (BoutSpikeTimes < GapOffset))) - INOnset;            
        end
    end
end

% Now to calculate firing rates in a gap window equal to the shortest last
% gap
INIndex = 0;
MinLastGap = min(Results(find(Position(:,1) == -1),4));

for i = 1:length(Neural_INR.NoofINs),
    if (Neural_INR.NoofINs(i) > 0)
        BoutSpikeTimes = Neural_INR.BoutDetails(i).SpikeTimes;
        INs = Neural_INR.INs{i};
        for j = 1:length(INs),
            INIndex = INIndex + 1;
            Position(INIndex, :) = [(j - length(INs) - 1) j length(INs) 1];
            GapOnset = Neural_INR.BoutDetails(i).offsets(INs(j)) - PreMotorLatency;
            GapOffset = GapOnset + MinLastGap;
            
            Results(INIndex, 7) = length(find((BoutSpikeTimes >= GapOnset) & (BoutSpikeTimes < GapOffset)))/(GapOffset - GapOnset);
            Results(INIndex, 8) = GapOffset - GapOnset;
            
            TempPST = histc(BoutSpikeTimes, (GapOnset:BinSize:GapOffset))/BinSize;
            GapPST(INIndex,:) = conv(TempPST(1:end-1), GaussWin, 'same');
            
            SpikeTimes{4}{INIndex} = BoutSpikeTimes(find((BoutSpikeTimes >= GapOnset) & (BoutSpikeTimes < GapOffset))) - GapOnset;
        end
    end
end

for i = 1:size(Neural_INR.WithinBoutNoofINs, 1),
    if (Neural_INR.WithinBoutNoofINs(i,1) > 0)
        BoutSpikeTimes = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes;
        INs = Neural_INR.WithinBoutINs{i};
        for j = 1:length(INs),
            INIndex = INIndex + 1;
            Position(INIndex, :) = [(j - length(INs) - 1) j length(INs) 0];
            GapOnset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(INs(j)) - PreMotorLatency;
            GapOffset = GapOnset + MinLastGap;
            
            Results(INIndex, 7) = length(find((BoutSpikeTimes >= GapOnset) & (BoutSpikeTimes < GapOffset)))/(GapOffset - GapOnset);
            Results(INIndex, 8) = GapOffset - GapOnset;
            
            TempPST = histc(BoutSpikeTimes, (GapOnset:BinSize:GapOffset))/BinSize;
            GapPST(INIndex,:) = conv(TempPST(1:end-1), GaussWin, 'same');
            
            SpikeTimes{4}{INIndex} = BoutSpikeTimes(find((BoutSpikeTimes >= GapOnset) & (BoutSpikeTimes < GapOffset))) - GapOnset;
        end
    end
end

% Now to warp spike times to median durations
MedianDurs = median(Results(:, [2 4 6]));

for i = 1:length(SpikeTimes{1}),
    WarpedSpikeTimes{1}{i} = SpikeTimes{1}{i} * MedianDurs(1)/Results(i,2);
    TempPST = histc(WarpedSpikeTimes{1}{i}, (0:BinSize:MedianDurs(1)))/BinSize;
    WarpedPST{1}(i,:) = conv(TempPST(1:end-1), GaussWin, 'same');
    
    WarpedSpikeTimes{2}{i} = SpikeTimes{2}{i} * MedianDurs(2)/Results(i,4);
    TempPST = histc(WarpedSpikeTimes{2}{i}, (0:BinSize:MedianDurs(2)))/BinSize;
    WarpedPST{2}(i,:) = conv(TempPST(1:end-1), GaussWin, 'same');
    
    WarpedSpikeTimes{3}{i} = SpikeTimes{3}{i} * MedianDurs(3)/Results(i,6);
    TempPST = histc(WarpedSpikeTimes{3}{i}, (0:BinSize:MedianDurs(3)))/BinSize;
    WarpedPST{3}(i,:) = conv(TempPST(1:end-1), GaussWin, 'same');
end

% Now for the stats on firing rate correlations
[r, p] = corrcoef(Results(:,1), Position(:,1));
Stats.INFiringRateCorr_r = r(1,2);
Stats.INFiringRateCorr_p = p(1,2);

[r, p] = corrcoef(Results(:,3), Position(:,1));
Stats.GapFiringRateCorr_r = r(1,2);
Stats.GapFiringRateCorr_p = p(1,2);

[r, p] = corrcoef(Results(:,5), Position(:,1));
Stats.INGapFiringRateCorr_r = r(1,2);
Stats.INGapFiringRateCorr_p = p(1,2);

[r, p] = corrcoef(Results(:,7), Position(:,1));
Stats.MinGapFiringRateCorr_r = r(1,2);
Stats.MinGapFiringRateCorr_p = p(1,2);

MinDurs = min(Results(:, [2 4 6]));
FullPST = PST;
clear PST;
for i = 1:3,
    for j = 1:length(FullPST{i}),
        PST{i}(j,:) = FullPST{i}{j}(1:floor(MinDurs(i)/BinSize));
    end
end

% % Now for the stats on firing pattern similarity - min duration
% for i = 1:3,
%    for j = 1:size(PST{i},1),
%         FirstINIndices = find((Position(:,2) == 1));
%         if (Position(j,2) == 1)
%             FirstINIndices(find(FirstINIndices == j)) = [];
%         end
%         LastINIndices = find((Position(:,1) == ComparisonIN));
%         if ((j-i-1) == ComparisonIN)
%             LastINIndices(find(LastINIndices == j)) = [];
%         end
%         
%         LastPST = mean(PST{i}(LastINIndices, :));
%         FirstPST = mean(PST{i}(FirstINIndices, :));
%         
%         TempCorr = corrcoef((LastPST - mean(LastPST)), (PST{i}(j,:) - mean(PST{i}(j,:))));
%         LastINCorr{i}(j) = TempCorr(1,2);
% 
%         TempCorr = corrcoef((FirstPST - mean(FirstPST)), (PST{i}(j,:) - mean(PST{i}(j,:))));
%         FirstINCorr{i}(j) = TempCorr(1,2);
%    end
% end
% 
% % Now for the stats on firing pattern similarity - min duration
% for i = 1:3,
%    for j = 1:size(WarpedPST{i},1),
%         FirstINIndices = find((Position(:,2) == 1));
%         if (Position(j,2) == 1)
%             FirstINIndices(find(FirstINIndices == j)) = [];
%         end
%         LastINIndices = find((Position(:,1) == ComparisonIN));
%         if ((j-i-1) == ComparisonIN)
%             LastINIndices(find(LastINIndices == j)) = [];
%         end
%         
%         WLastPST = mean(WarpedPST{i}(LastINIndices, :));
%         WFirstPST = mean(WarpedPST{i}(FirstINIndices, :));
%         
%         TempCorr = corrcoef((WLastPST - mean(WLastPST)), (WarpedPST{i}(j,:) - mean(WarpedPST{i}(j,:))));
%         WLastINCorr{i}(j) = TempCorr(1,2);
% 
%         TempCorr = corrcoef((WFirstPST - mean(WFirstPST)), (WarpedPST{i}(j,:) - mean(WarpedPST{i}(j,:))));
%         WFirstINCorr{i}(j) = TempCorr(1,2);
%    end
% end

disp('Finished');