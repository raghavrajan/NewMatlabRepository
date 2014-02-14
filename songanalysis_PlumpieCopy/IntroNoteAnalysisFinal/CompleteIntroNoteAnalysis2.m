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

MinDurs = min(Results(:, [2 4 6]));
FullPST = PST;
clear PST;
for i = 1:3,
    for j = 1:length(FullPST{i}),
        PST{i}(j,:) = FullPST{i}{j}(1:floor(MinDurs(i)/BinSize));
    end
end


% Now for the stats on firing pattern similarity - min duration
for i = min(Position(:,3)):1:max(Position(:,3)),
    Indices = find((Position(:,3) == i) & (Position(:,1) == -1));
    TrialNos(i) = length(Indices);
    if (length(Indices) >= 2)
        for j = 1:i,
            if (j == 1)
                FirstINIndices = find((Position(:,2) == 1) & (Position(:,3) ~=i));
            else
                FirstINIndices = find((Position(:,2) == 1));
            end
            clear FirstINPST LastINPST SpecificINPST RandomINPST;
            FirstINPST = mean(PST{1}(FirstINIndices, :));
            FirstGapPST = mean(PST{2}(FirstINIndices, :));
            FirstINGapPST = mean(PST{3}(FirstINIndices, :));
            
            if ((j-i-1) == ComparisonIN)
                LastINIndices = find((Position(:,1) == ComparisonIN) & (Position(:,3) ~=i));
            else
                LastINIndices = find((Position(:,1) == ComparisonIN));
            end
            
            if (length(LastINIndices) > 1)
                LastINPST = mean(PST{1}(LastINIndices, :));
                LastGapPST = mean(PST{2}(LastINIndices, :));
                LastINGapPST = mean(PST{3}(LastINIndices, :));
            else
                LastINPST = (PST{1}(LastINIndices, :));
                LastGapPST = (PST{2}(LastINIndices, :));
                LastINGapPST = (PST{3}(LastINIndices, :));
            end

            SpecificINIndices = find((Position(:,3) == i) & (Position(:,1) == (j - i - 1)));
            SpecificINPST = mean(PST{1}(SpecificINIndices, :));
            SpecificGapPST = mean(PST{2}(SpecificINIndices, :));
            SpecificINGapPST = mean(PST{3}(SpecificINIndices, :));

            TempCorr = corrcoef((LastINPST - mean(LastINPST)), (SpecificINPST - mean(SpecificINPST)));
            LastINCorr{i}(j) = TempCorr(1,2);

            TempCorr = corrcoef((LastGapPST - mean(LastGapPST)), (SpecificGapPST - mean(SpecificGapPST)));
            LastGapCorr{i}(j) = TempCorr(1,2);
            
            TempCorr = corrcoef((LastINGapPST - mean(LastINGapPST)), (SpecificINGapPST - mean(SpecificINGapPST)));
            LastINGapCorr{i}(j) = TempCorr(1,2);
            
            TempCorr = corrcoef((FirstINPST - mean(FirstINPST)), (SpecificINPST - mean(SpecificINPST)));
            FirstINCorr{i}(j) = TempCorr(1,2);
            
            TempCorr = corrcoef((FirstGapPST - mean(FirstGapPST)), (SpecificGapPST - mean(SpecificGapPST)));
            FirstGapCorr{i}(j) = TempCorr(1,2);
            
            TempCorr = corrcoef((FirstINGapPST - mean(FirstINGapPST)), (SpecificINGapPST - mean(SpecificINGapPST)));
            FirstINGapCorr{i}(j) = TempCorr(1,2);
        end
        for j = 1:1000,
            Indices = find(Position(:,3) == i);
            RandomINIndices = randperm(length(Indices));
            RandomINIndices = Indices(RandomINIndices(1:length(Indices)/i));
            RandomINPST = mean(PST{1}(RandomINIndices, :));
            RandomGapPST = mean(PST{2}(RandomINIndices, :));
            RandomINGapPST = mean(PST{3}(RandomINIndices, :));

            TempCorr = corrcoef((LastINPST - mean(LastINPST)), (RandomINPST - mean(RandomINPST)));
            RandomLastINCorr{i}(j) = TempCorr(1,2);

            TempCorr = corrcoef((LastGapPST - mean(LastGapPST)), (RandomGapPST - mean(RandomGapPST)));
            RandomLastGapCorr{i}(j) = TempCorr(1,2);

            TempCorr = corrcoef((LastINGapPST - mean(LastINGapPST)), (RandomINGapPST - mean(RandomINGapPST)));
            RandomLastINGapCorr{i}(j) = TempCorr(1,2);

            TempCorr = corrcoef((FirstINPST - mean(FirstINPST)), (RandomINPST - mean(RandomINPST)));
            RandomFirstINCorr{i}(j) = TempCorr(1,2);

            TempCorr = corrcoef((FirstGapPST - mean(FirstGapPST)), (RandomGapPST - mean(RandomGapPST)));
            RandomFirstGapCorr{i}(j) = TempCorr(1,2);

            TempCorr = corrcoef((FirstINGapPST - mean(FirstINGapPST)), (RandomINGapPST - mean(RandomINGapPST)));
            RandomFirstINGapCorr{i}(j) = TempCorr(1,2);
        end
        RandomLastINCorr{i} = sort(RandomLastINCorr{i});
        RandomLastINCorr{i} = [RandomLastINCorr{i}(round(0.05*1000)) RandomLastINCorr{i}(round(0.95*1000))];
        
        RandomLastGapCorr{i} = sort(RandomLastGapCorr{i});
        RandomLastGapCorr{i} = [RandomLastGapCorr{i}(round(0.05*1000)) RandomLastGapCorr{i}(round(0.95*1000))];
        
        RandomLastINGapCorr{i} = sort(RandomLastINGapCorr{i});
        RandomLastINGapCorr{i} = [RandomLastINGapCorr{i}(round(0.05*1000)) RandomLastINGapCorr{i}(round(0.95*1000))];
        
        RandomFirstINCorr{i} = sort(RandomFirstINCorr{i});
        RandomFirstINCorr{i} = [RandomFirstINCorr{i}(round(0.05*1000)) RandomFirstINCorr{i}(round(0.95*1000))];
        
        RandomFirstGapCorr{i} = sort(RandomFirstGapCorr{i});
        RandomFirstGapCorr{i} = [RandomFirstGapCorr{i}(round(0.05*1000)) RandomFirstGapCorr{i}(round(0.95*1000))];
        
        RandomFirstINGapCorr{i} = sort(RandomFirstINGapCorr{i});
        RandomFirstINGapCorr{i} = [RandomFirstINGapCorr{i}(round(0.05*1000)) RandomFirstINGapCorr{i}(round(0.95*1000))];
    end
end

% Now for the stats on firing pattern similarity - warped psts

for i = min(Position(:,3)):1:max(Position(:,3)),
    Indices = find((Position(:,3) == i) & (Position(:,1) == -1));
    if (length(Indices) >= 2)
        for j = 1:i,
            if (j == 1)
                FirstINIndices = find((Position(:,2) == 1) & (Position(:,3) ~=i));
            else
                FirstINIndices = find((Position(:,2) == 1));
            end
            clear FirstINPST LastINPST SpecificINPST RandomINPST;
            FirstINPST = mean(WarpedPST{1}(FirstINIndices, :));
            FirstGapPST = mean(WarpedPST{2}(FirstINIndices, :));
            FirstINGapPST = mean(WarpedPST{3}(FirstINIndices, :));
            
            if ((j-i-1) == ComparisonIN)
                LastINIndices = find((Position(:,1) == ComparisonIN) & (Position(:,3) ~=i));
            else
                LastINIndices = find((Position(:,1) == ComparisonIN));
            end
            
            if (length(LastINIndices) > 1)
                LastINPST = mean(WarpedPST{1}(LastINIndices, :));
                LastGapPST = mean(WarpedPST{2}(LastINIndices, :));
                LastINGapPST = mean(WarpedPST{3}(LastINIndices, :));
            else
                LastINPST = (WarpedPST{1}(LastINIndices, :));
                LastGapPST = (WarpedPST{2}(LastINIndices, :));
                LastINGapPST = (WarpedPST{3}(LastINIndices, :));
            end

            SpecificINIndices = find((Position(:,3) == i) & (Position(:,1) == (j - i - 1)));
            SpecificINPST = mean(WarpedPST{1}(SpecificINIndices, :));
            SpecificGapPST = mean(WarpedPST{2}(SpecificINIndices, :));
            SpecificINGapPST = mean(WarpedPST{3}(SpecificINIndices, :));

            TempCorr = corrcoef((LastINPST - mean(LastINPST)), (SpecificINPST - mean(SpecificINPST)));
            WLastINCorr{i}(j) = TempCorr(1,2);

            TempCorr = corrcoef((LastGapPST - mean(LastGapPST)), (SpecificGapPST - mean(SpecificGapPST)));
            WLastGapCorr{i}(j) = TempCorr(1,2);
            
            TempCorr = corrcoef((LastINGapPST - mean(LastINGapPST)), (SpecificINGapPST - mean(SpecificINGapPST)));
            WLastINGapCorr{i}(j) = TempCorr(1,2);
            
            TempCorr = corrcoef((FirstINPST - mean(FirstINPST)), (SpecificINPST - mean(SpecificINPST)));
            WFirstINCorr{i}(j) = TempCorr(1,2);
            
            TempCorr = corrcoef((FirstGapPST - mean(FirstGapPST)), (SpecificGapPST - mean(SpecificGapPST)));
            WFirstGapCorr{i}(j) = TempCorr(1,2);
            
            TempCorr = corrcoef((FirstINGapPST - mean(FirstINGapPST)), (SpecificINGapPST - mean(SpecificINGapPST)));
            WFirstINGapCorr{i}(j) = TempCorr(1,2);
        end
        for j = 1:1000,
            Indices = find(Position(:,3) == i);
            RandomINIndices = randperm(length(Indices));
            RandomINIndices = Indices(RandomINIndices(1:length(Indices)/i));
            RandomINPST = mean(WarpedPST{1}(RandomINIndices, :));
            RandomGapPST = mean(WarpedPST{2}(RandomINIndices, :));
            RandomINGapPST = mean(WarpedPST{3}(RandomINIndices, :));

            TempCorr = corrcoef((LastINPST - mean(LastINPST)), (RandomINPST - mean(RandomINPST)));
            WRandomLastINCorr{i}(j) = TempCorr(1,2);

            TempCorr = corrcoef((LastGapPST - mean(LastGapPST)), (RandomGapPST - mean(RandomGapPST)));
            WRandomLastGapCorr{i}(j) = TempCorr(1,2);

            TempCorr = corrcoef((LastINGapPST - mean(LastINGapPST)), (RandomINGapPST - mean(RandomINGapPST)));
            WRandomLastINGapCorr{i}(j) = TempCorr(1,2);

            TempCorr = corrcoef((FirstINPST - mean(FirstINPST)), (RandomINPST - mean(RandomINPST)));
            WRandomFirstINCorr{i}(j) = TempCorr(1,2);

            TempCorr = corrcoef((FirstGapPST - mean(FirstGapPST)), (RandomGapPST - mean(RandomGapPST)));
            WRandomFirstGapCorr{i}(j) = TempCorr(1,2);

            TempCorr = corrcoef((FirstINGapPST - mean(FirstINGapPST)), (RandomINGapPST - mean(RandomINGapPST)));
            WRandomFirstINGapCorr{i}(j) = TempCorr(1,2);
        end
        WRandomLastINCorr{i} = sort(WRandomLastINCorr{i});
        WRandomLastINCorr{i} = [WRandomLastINCorr{i}(round(0.05*1000)) WRandomLastINCorr{i}(round(0.95*1000))];
        
        WRandomLastGapCorr{i} = sort(WRandomLastGapCorr{i});
        WRandomLastGapCorr{i} = [WRandomLastGapCorr{i}(round(0.05*1000)) WRandomLastGapCorr{i}(round(0.95*1000))];
        
        WRandomLastINGapCorr{i} = sort(WRandomLastINGapCorr{i});
        WRandomLastINGapCorr{i} = [WRandomLastINGapCorr{i}(round(0.05*1000)) WRandomLastINGapCorr{i}(round(0.95*1000))];
        
        WRandomFirstINCorr{i} = sort(WRandomFirstINCorr{i});
        WRandomFirstINCorr{i} = [WRandomFirstINCorr{i}(round(0.05*1000)) WRandomFirstINCorr{i}(round(0.95*1000))];
        
        WRandomFirstGapCorr{i} = sort(WRandomFirstGapCorr{i});
        WRandomFirstGapCorr{i} = [WRandomFirstGapCorr{i}(round(0.05*1000)) WRandomFirstGapCorr{i}(round(0.95*1000))];
        
        WRandomFirstINGapCorr{i} = sort(WRandomFirstINGapCorr{i});
        WRandomFirstINGapCorr{i} = [WRandomFirstINGapCorr{i}(round(0.05*1000)) WRandomFirstINGapCorr{i}(round(0.95*1000))];
    end
end

MaxPlots = max(Position(:,3));

ValidTrials = find(TrialNos >= 2);

AllFirstINCorr = ones(length(ValidTrials), MaxPlots)*NaN;
AllFirstGapCorr = ones(length(ValidTrials), MaxPlots)*NaN;
AllFirstINGapCorr = ones(length(ValidTrials), MaxPlots)*NaN;

AllLastINCorr = ones(length(ValidTrials), MaxPlots)*NaN;
AllLastGapCorr = ones(length(ValidTrials), MaxPlots)*NaN;
AllLastINGapCorr = ones(length(ValidTrials), MaxPlots)*NaN;

AllRandomINCorr = ones(length(ValidTrials), MaxPlots)*NaN;
AllRandomGapCorr = ones(length(ValidTrials), MaxPlots)*NaN;
AllRandomINGapCorr = ones(length(ValidTrials), MaxPlots)*NaN;

for j = 1:length(ValidTrials),
    for k = 1:ValidTrials(j),
        AllFirstINCorr(j, (MaxPlots - (ValidTrials(j) - k))) = FirstINCorr{ValidTrials(j)}(k);
        AllLastINCorr(j, (MaxPlots - (ValidTrials(j) - k))) = LastINCorr{ValidTrials(j)}(k);
        AllFFWinValue(j, (MaxPlots - (ValidTrials(j) - k))) = INGapFFWindowValue{ValidTrials(j)}(k);
    end
end

for j = 1:size(AllFR,2),
    TempAllFR(j,:) = [mean(AllFR(find(~isnan(AllFR(:,j))),j)) std(AllFR(find(~isnan(AllFR(:,j))),j))/sqrt(length(find(~isnan(AllFR(:,j))))) j];
    TempAllFirstINCorr(j,:) = [mean(AllFirstINCorr(find(~isnan(AllFirstINCorr(:,j))),j)) std(AllFirstINCorr(find(~isnan(AllFirstINCorr(:,j))),j))/sqrt(length(find(~isnan(AllFirstINCorr(:,j)))))];
    TempAllLastINCorr(j,:) = [mean(AllLastINCorr(find(~isnan(AllLastINCorr(:,j))),j)) std(AllLastINCorr(find(~isnan(AllLastINCorr(:,j))),j))/sqrt(length(find(~isnan(AllLastINCorr(:,j)))))];
    TempAllFFWinValue(j,:) = [mean(AllFFWinValue(find(~isnan(AllFFWinValue(:,j))),j)) std(AllFFWinValue(find(~isnan(AllFFWinValue(:,j))),j))/sqrt(length(find(~isnan(AllFFWinValue(:,j)))))];
    TempAllINFR(j,:) = [mean(AllINFR(find(~isnan(AllINFR(:,j))),j)) std(AllINFR(find(~isnan(AllINFR(:,j))),j))/sqrt(length(find(~isnan(AllINFR(:,j))))) j];
    TempAllGapFR(j,:) = [mean(AllGapFR(find(~isnan(AllGapFR(:,j))),j)) std(AllGapFR(find(~isnan(AllGapFR(:,j))),j))/sqrt(length(find(~isnan(AllGapFR(:,j))))) j];
end

disp('Finished');