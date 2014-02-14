function [] = IntroNoteWarpedRasters(Neural_INR, BinSize, PreMotorLatency, MinNumber, MinPlotNumber)

FFWindow = 0.03;

Width = 0.005;
GaussianLen = 2;
IFRFs = 1/BinSize;
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (IFRFs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * IFRFs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * IFRFs) * (Width * IFRFs)));

PreTime = 0.2;

[INDetails] = IntroNoteAnalysisGetSyllGapSpikeTimes(Neural_INR{1}, 0, BinSize);

PosIndex = 0;
for i = max([INDetails.INs.Position]):-1:min([INDetails.INs.Position]),
    Trials = find([INDetails.INs.Position] == i);
    if (~isempty(Trials))
        PosIndex = PosIndex + 1;
        Position(PosIndex).Index = i;
        Position(PosIndex).Trials = Trials;
        Position(PosIndex).MedianSyllDur = median([INDetails.INs(Position(PosIndex).Trials).INDur]);
        Position(PosIndex).MedianGapDur = median([INDetails.INs(Position(PosIndex).Trials).GapDur]);
        Position(PosIndex).MedianSyllGapDur = median([INDetails.INs(Position(PosIndex).Trials).INGapDur]);
    end
end

MedianSyllDurs = [Position.MedianSyllDur];
MedianGapDurs = [Position.MedianGapDur];
MedianSyllDurs = ones(size(MedianSyllDurs))*MedianSyllDurs(1);
MedianGapDurs = ones(size(MedianGapDurs))*MedianGapDurs(1);

MedianFirstSyllDur = median([INDetails.FirstSylls.FirstSyllDur]);


for i = 1:max(abs([Position.Index])),

    Edges = MedianFirstSyllDur:-BinSize:-(PreTime + sum(MedianSyllDurs(1:i)) + sum(MedianGapDurs(1:i)));
    Edges = fliplr(Edges);
    FFEdges = (-(PreTime + sum(MedianSyllDurs(1:i)) + sum(MedianGapDurs(1:i))) + FFWindow/2):0.001:(MedianFirstSyllDur - FFWindow/2);
    RasterPST(i).NoofINs = i; 
    RasterPST(i).Raster = [];
    RasterPST(i).PST = [];
    RasterPST(i).Edges = Edges;
    RasterPST(i).FFEdges = FFEdges;

    TrialIndex = 0;
    
    % Find the number of trials with 'i' number of INs - both from the
    % beginning of bouts and within bouts
    
    Trials = length(find((Neural_INR{1}.NoofINs == i))) + length(find(Neural_INR{1}.WithinBoutNoofINs(:,1) == i));
    
    if (Trials > MinNumber)
        
        % Now first get all the data associated with the beginning of bouts
        
        BegTrials = find(Neural_INR{1}.NoofINs == i);
        for j = 1:length(BegTrials),
            TrialSpikeTimes = [];
            TrialIndex = TrialIndex + 1;
            BoutSpikeTimes = Neural_INR{1}.BoutDetails(BegTrials(j)).SpikeTimes + PreMotorLatency;
            BoutSpikeTimes = BoutSpikeTimes(:);
            INs = Neural_INR{1}.INs{BegTrials(j)};
            % First get pre-time spikes
            PreTimeEnd = -(sum(MedianSyllDurs(1:i)) + sum(MedianGapDurs(1:i)));
            PreTimeStart = PreTimeEnd - PreTime;
            
            % First get the spikes occuring in the short pre-time period
            % before the first intro note
            TempSpikeIndices = find((BoutSpikeTimes >= (Neural_INR{1}.BoutDetails(BegTrials(j)).onsets(INs(1)) - PreTime)) & (BoutSpikeTimes < Neural_INR{1}.BoutDetails(BegTrials(j)).onsets(INs(1))));
            TempSpikeIndices = TempSpikeIndices(:);
            
            if (~isempty(TempSpikeIndices))
                TrialSpikeTimes = [TrialSpikeTimes; (BoutSpikeTimes(TempSpikeIndices) - (Neural_INR{1}.BoutDetails(BegTrials(j)).onsets(INs(1)) - PreTime) + PreTimeStart)];
            end
            
            % Now get the spikes during each IN and during each gap and
            % warp it to the median IN length and median gap length -
            % starting at the end, so that the end is the last IN
            
            for k = length(INs):-1:1,
                % First get spikes associated with the gap
                NextSyllONTime = Neural_INR{1}.BoutDetails(BegTrials(j)).onsets(INs(k) + 1);
                SyllOFFTime = Neural_INR{1}.BoutDetails(BegTrials(j)).offsets(INs(k));
                TempSpikeIndices = find((BoutSpikeTimes >= SyllOFFTime) & (BoutSpikeTimes < NextSyllONTime));
                TempSpikeIndices = TempSpikeIndices(:);
                if (~isempty(TempSpikeIndices))
                    WarpedSpikeTimes = (BoutSpikeTimes(TempSpikeIndices) - SyllOFFTime) * MedianGapDurs(length(INs) - k + 1) / (NextSyllONTime - SyllOFFTime);
                
                    TrialSpikeTimes = [TrialSpikeTimes; (WarpedSpikeTimes - sum(MedianGapDurs(1:(length(INs) - k + 1))) - sum(MedianSyllDurs(1:(length(INs) - k))))];
                end
                
                % Now get the spikes associated with the IN
                SyllONTime = Neural_INR{1}.BoutDetails(BegTrials(j)).onsets(INs(k));
                SyllOFFTime = Neural_INR{1}.BoutDetails(BegTrials(j)).offsets(INs(k));
                TempSpikeIndices = find((BoutSpikeTimes >= SyllONTime) & (BoutSpikeTimes < SyllOFFTime));
                TempSpikeIndices = TempSpikeIndices(:);
                if (~isempty(TempSpikeIndices))
                    WarpedSpikeTimes = (BoutSpikeTimes(TempSpikeIndices) - SyllONTime) * MedianSyllDurs(length(INs) - k + 1) / (SyllOFFTime - SyllONTime);
                    TrialSpikeTimes = [TrialSpikeTimes; (WarpedSpikeTimes - sum(MedianGapDurs(1:(length(INs) - k + 1))) - sum(MedianSyllDurs(1:(length(INs) - k + 1))))];
                end

            end
            % Now get the spikes associated with the first syllable
            SyllONTime = Neural_INR{1}.BoutDetails(BegTrials(j)).onsets(INs(end) + 1);
            SyllOFFTime = Neural_INR{1}.BoutDetails(BegTrials(j)).offsets(INs(end) + 1);
            TempSpikeIndices = find((BoutSpikeTimes >= SyllONTime) & (BoutSpikeTimes < SyllOFFTime));
            TempSpikeIndices = TempSpikeIndices(:);
            if (~isempty(TempSpikeIndices))
                WarpedSpikeTimes = (BoutSpikeTimes(TempSpikeIndices) - SyllONTime) * MedianFirstSyllDur / (SyllOFFTime - SyllONTime);
                TrialSpikeTimes = [TrialSpikeTimes; WarpedSpikeTimes];
            end
            
            TrialSpikeTimes = sort(TrialSpikeTimes);
            
            RasterPST(i).TrialSpikeTimes{TrialIndex} = TrialSpikeTimes;
            RasterPST(i).Raster = [RasterPST(i).Raster; [TrialSpikeTimes ones(size(TrialSpikeTimes))*TrialIndex]];
            if (~isempty(TrialSpikeTimes))
                TempTrialPST = histc(TrialSpikeTimes, Edges)/BinSize;
                TempTrialPST = TempTrialPST(:)';
                RasterPST(i).PST(TrialIndex, :) = TempTrialPST;
            else
                RasterPST(i).PST(TrialIndex, :) = zeros(size(Edges));
            end
            for FFIndex = 1:length(FFEdges),
                RasterPST(i).FF(TrialIndex, FFIndex) = length(find((TrialSpikeTimes >= (FFEdges(FFIndex) - FFWindow/2)) & (TrialSpikeTimes < (FFEdges(FFIndex) + FFWindow/2))));
            end
        end
        
        % Now get the data associated with within bouts
        
        WithinBoutTrials = find(Neural_INR{1}.WithinBoutNoofINs(:,1) == i);
        for j = 1:length(WithinBoutTrials),
            TrialSpikeTimes = [];
            
            TrialIndex = TrialIndex + 1;
            BoutIndex = Neural_INR{1}.WithinBoutINBoutIndices(WithinBoutTrials(j));
            BoutSpikeTimes = Neural_INR{1}.BoutDetails(BoutIndex).SpikeTimes + PreMotorLatency;
            BoutSpikeTimes = BoutSpikeTimes(:);
            INs = Neural_INR{1}.WithinBoutINs{WithinBoutTrials(j)};
            
            % First get pre-time spikes
            PreTimeEnd = -(sum(MedianSyllDurs(1:i)) + sum(MedianGapDurs(1:i)));
            PreTimeStart = PreTimeEnd - PreTime;
            
            % First get the spikes occuring in the short pre-time period
            % before the first intro note
            TempSpikeIndices = find((BoutSpikeTimes >= (Neural_INR{1}.BoutDetails(BoutIndex).onsets(INs(1)) - PreTime)) & (BoutSpikeTimes < Neural_INR{1}.BoutDetails(BoutIndex).onsets(INs(1))));
            TempSpikeIndices = TempSpikeIndices(:);
            if (~isempty(TempSpikeIndices))
                TrialSpikeTimes = [TrialSpikeTimes; (BoutSpikeTimes(TempSpikeIndices) - (Neural_INR{1}.BoutDetails(BoutIndex).onsets(INs(1)) - PreTime) + PreTimeStart)];
            end
            
            % Now get the spikes during each IN and during each gap and
            % warp it to the median IN length and median gap length -
            % starting at the end, so that the end is the last IN
            
            for k = length(INs):-1:1,
                % First get spikes associated with the gap
                NextSyllONTime = Neural_INR{1}.BoutDetails(BoutIndex).onsets(INs(k) + 1);
                SyllOFFTime = Neural_INR{1}.BoutDetails(BoutIndex).offsets(INs(k));
                TempSpikeIndices = find((BoutSpikeTimes >= SyllOFFTime) & (BoutSpikeTimes < NextSyllONTime));
                TempSpikeIndices = TempSpikeIndices(:);
                if (~isempty(TempSpikeIndices))
                    WarpedSpikeTimes = (BoutSpikeTimes(TempSpikeIndices) - SyllOFFTime) * MedianGapDurs(length(INs) - k + 1) / (NextSyllONTime - SyllOFFTime);
                    TrialSpikeTimes = [TrialSpikeTimes; (WarpedSpikeTimes - sum(MedianGapDurs(1:(length(INs) - k + 1))) - sum(MedianSyllDurs(1:(length(INs) - k))))];
                end
                
                % Now get the spikes associated with the IN
                SyllONTime = Neural_INR{1}.BoutDetails(BoutIndex).onsets(INs(k));
                SyllOFFTime = Neural_INR{1}.BoutDetails(BoutIndex).offsets(INs(k));
                TempSpikeIndices = find((BoutSpikeTimes >= SyllONTime) & (BoutSpikeTimes < SyllOFFTime));
                TempSpikeIndices = TempSpikeIndices(:);
                if (~isempty(TempSpikeIndices))
                    WarpedSpikeTimes = (BoutSpikeTimes(TempSpikeIndices) - SyllONTime) * MedianSyllDurs(length(INs) - k + 1) / (SyllOFFTime - SyllONTime);
                    TrialSpikeTimes = [TrialSpikeTimes; (WarpedSpikeTimes - sum(MedianGapDurs(1:(length(INs) - k + 1))) - sum(MedianSyllDurs(1:(length(INs) - k + 1))))];
                end
            end
            
            % Now get the spikes associated with the first syllable
            SyllONTime = Neural_INR{1}.BoutDetails(BoutIndex).onsets(INs(end) + 1);
            SyllOFFTime = Neural_INR{1}.BoutDetails(BoutIndex).offsets(INs(end) + 1);
            TempSpikeIndices = find((BoutSpikeTimes >= SyllONTime) & (BoutSpikeTimes < SyllOFFTime));
            TempSpikeIndices = TempSpikeIndices(:);
            if (~isempty(TempSpikeIndices))
                WarpedSpikeTimes = (BoutSpikeTimes(TempSpikeIndices) - SyllONTime) * MedianFirstSyllDur / (SyllOFFTime - SyllONTime);
                TrialSpikeTimes = [TrialSpikeTimes; WarpedSpikeTimes];
            end
            
            TrialSpikeTimes = sort(TrialSpikeTimes);
            
            RasterPST(i).TrialSpikeTimes{TrialIndex} = TrialSpikeTimes;
            RasterPST(i).Raster = [RasterPST(i).Raster; [TrialSpikeTimes ones(size(TrialSpikeTimes))*TrialIndex]];
            if (~isempty(TrialSpikeTimes))
                RasterPST(i).PST(TrialIndex, :) = histc(TrialSpikeTimes, Edges)/BinSize;
            else
                RasterPST(i).PST(TrialIndex, :) = zeros(size(Edges));
            end
            for FFIndex = 1:length(FFEdges),
                RasterPST(i).FF(TrialIndex, FFIndex) = length(find((TrialSpikeTimes >= (FFEdges(FFIndex) - FFWindow/2)) & (TrialSpikeTimes < (FFEdges(FFIndex) + FFWindow/2))));
            end
        end
    end
end

Colours = 'rbk';
figure;
RasterIncrement = 0;
ValidNums = [];
for i = 1:length(RasterPST),
    if (size(RasterPST(i).PST, 1) >= MinPlotNumber)
        ValidNums = [ValidNums; i];
        subplot(8,1,[1:4]);
        PlotRaster(RasterPST(i).Raster, Colours(mod(i-1,length(Colours)) + 1), 0.25, RasterIncrement, MinPlotNumber);
        hold on;
        RasterIncrement = RasterIncrement + MinPlotNumber + 1;
        
        subplot(8,1,[5:6]);
        hold on;
        plot(RasterPST(i).Edges, conv(mean(RasterPST(i).PST), GaussWin, 'same'), Colours(mod(i-1, length(Colours)) + 1));
        
        subplot(8,1,[7:8]);
        hold on;
        plot(RasterPST(i).FFEdges, var(RasterPST(i).FF)./mean(RasterPST(i).FF), Colours(mod(i-1, length(Colours)) + 1));
    end
end
subplot(8,1,[1:4]);
axis tight;
temp = axis;
axis([RasterPST(ValidNums(end)).Edges(1) RasterPST(ValidNums(end)).Edges(end) 0 1.05*temp(4)]);

for i = 1:ValidNums(end),
    Temp = find(ValidNums >= i, 1, 'first');
    YValues = [(Temp-1)*(round(1+MinPlotNumber))+0.5 (Temp-1)*(1+MinPlotNumber)+0.5 temp(4) temp(4)];
    fill([-(sum(MedianSyllDurs(1:i)) + sum(MedianGapDurs(1:i))) -(sum(MedianGapDurs(1:i)) + sum(MedianSyllDurs(1:i-1))) fliplr([-(sum(MedianSyllDurs(1:i)) + sum(MedianGapDurs(1:i))) -(sum(MedianGapDurs(1:i)) + sum(MedianSyllDurs(1:i-1)))])], YValues, 'k', 'EdgeColor', 'none', 'FaceColor', [0.75 0.75 0.75], 'FaceAlpha', 0.5);
end
YValues = [0.5 0.5 temp(4) temp(4)];
fill([0 MedianFirstSyllDur MedianFirstSyllDur 0], YValues, 'k', 'EdgeColor', 'none', 'FaceColor', [0.75 0.75 0.75], 'FaceAlpha', 0.5);

subplot(8,1,[5:6]);
axis tight;
temp = axis;
axis([RasterPST(ValidNums(end)).Edges(1) RasterPST(ValidNums(end)).Edges(end) 0 1.05*temp(4)]);
for i = 1:ValidNums(end),
    fill([-(sum(MedianSyllDurs(1:i)) + sum(MedianGapDurs(1:i))) -(sum(MedianGapDurs(1:i)) + sum(MedianSyllDurs(1:i-1))) fliplr([-(sum(MedianSyllDurs(1:i)) + sum(MedianGapDurs(1:i))) -(sum(MedianGapDurs(1:i)) + sum(MedianSyllDurs(1:i-1)))])], [0.5 0.5 temp(4) temp(4)], 'k', 'EdgeColor', 'none', 'FaceColor', [0.75 0.75 0.75], 'FaceAlpha', 0.5);
end
YValues = [0.5 0.5 temp(4) temp(4)];
fill([0 MedianFirstSyllDur MedianFirstSyllDur 0], YValues, 'k', 'EdgeColor', 'none', 'FaceColor', [0.75 0.75 0.75], 'FaceAlpha', 0.5);

subplot(8,1,[7:8]);
axis tight;
temp = axis;
axis([RasterPST(ValidNums(end)).Edges(1) RasterPST(ValidNums(end)).Edges(end) 0 1.05*temp(4)]);
for i = 1:ValidNums(end),
    fill([-(sum(MedianSyllDurs(1:i)) + sum(MedianGapDurs(1:i))) -(sum(MedianGapDurs(1:i)) + sum(MedianSyllDurs(1:i-1))) fliplr([-(sum(MedianSyllDurs(1:i)) + sum(MedianGapDurs(1:i))) -(sum(MedianGapDurs(1:i)) + sum(MedianSyllDurs(1:i-1)))])], [0 0 temp(4) temp(4)], 'k', 'EdgeColor', 'none', 'FaceColor', [0.75 0.75 0.75], 'FaceAlpha', 0.5);
end
YValues = [0 0 temp(4) temp(4)];
fill([0 MedianFirstSyllDur MedianFirstSyllDur 0], YValues, 'k', 'EdgeColor', 'none', 'FaceColor', [0.75 0.75 0.75], 'FaceAlpha', 0.5);
hold on;
plot([RasterPST(ValidNums(end)).Edges(1) RasterPST(ValidNums(end)).Edges(end)], [1 1], 'k--');
disp('Finished');