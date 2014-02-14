function [PC] = IntroNoteUnWarpedRastersAllNeuronsBird(Neural_INR, BinSize, PreMotorLatency, MinNumber, MinPlotNumber, INNum)

PreTime = 0.075;
PostTime = 0.075;

for i = 1:length(Neural_INR),
    c;
    [INGapRaster{i}, INPST{i}, PST{i}, Position{i}] = IntroNoteUnWarpedRastersAllData(Neural_INR{i}, BinSize, MinNumber, MinPlotNumber, PreTime, PostTime);
end

TempPST = ones(INNum*(size(PST{i},2)-1), length(INGapRaster));

PSTSize = size(PST{i},2) - 1;

for i = 1:length(INGapRaster),
    if (~isempty(INGapRaster{i}{INNum}))
        for j = 1:INNum,
            TempPST(((j-1)*PSTSize + 1):j*PSTSize, i) = mean(INPST{i}{INNum}{j}(:,1:end-1));
        end
    else
        TempPST(:,i) = ones(PSTSize*INNum, 1) * NaN;
    end
end

if (nargin > 6)
    PC = varargin{1};
    Symbol = varargin{2};
end

Width = 0.005;
GaussianLen = 2;
IFRFs = 1/BinSize;
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (IFRFs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * IFRFs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * IFRFs) * (Width * IFRFs)));

PreTime = 0.2;

for i = 1:length(Neural_INR),
    [INDetails(i)] = IntroNoteAnalysisGetSyllGapSpikeTimes(Neural_INR{i}, 0, BinSize);
end

MinPos = 1000;
MaxPos = -1000;
for i = 1:length(INDetails),
    MaxPos = max([MaxPos max([INDetails(i).INs.Position])]);
    MinPos = min([MinPos min([INDetails(i).INs.Position])]);
end

PosIndex = 0;
for i = MaxPos:-1:MinPos,
    TempSyllDur = [];
    TempGapDur = [];
    TempSyllGapDur = [];
    for NeuronNo = 1:length(INDetails),
        Trials = find([INDetails(NeuronNo).INs.Position] == i);
        if (~isempty(Trials))
            TempSyllDur = [TempSyllDur [INDetails(NeuronNo).INs(Trials).INDur]];
            TempGapDur = [TempGapDur [INDetails(NeuronNo).INs(Trials).GapDur]];
            TempSyllGapDur = [TempSyllGapDur [INDetails(NeuronNo).INs(Trials).INGapDur]];
        end
    end
    if (~isempty(TempSyllDur))
        PosIndex = PosIndex + 1;
        Position(PosIndex).Index = i;
        Position(PosIndex).MedianSyllDur = median(TempSyllDur);
        Position(PosIndex).MedianGapDur = median(TempGapDur);
        Position(PosIndex).MedianSyllGapDur = median(TempSyllGapDur);
    end
end

MedianSyllDurs = [Position.MedianSyllDur];
MedianGapDurs = [Position.MedianGapDur];

for i = MaxPos:-1:MinPos,
    TempSyllDur = [];
    for NeuronNo = 1:length(INDetails),
        TempSyllDur = [TempSyllDur [INDetails(NeuronNo).FirstSylls.FirstSyllDur]];
    end
end
MedianFirstSyllDur = median(TempSyllDur);


for NeuronNo = 1:length(Neural_INR),
    for i = 1:max(abs([Position.Index])),

        Edges = MedianFirstSyllDur:-BinSize:-(PreTime + sum(MedianSyllDurs(1:i)) + sum(MedianGapDurs(1:i)));
        Edges = fliplr(Edges);

        NeuronData{NeuronNo}.RasterPST(i).NoofINs = i; 
        NeuronData{NeuronNo}.RasterPST(i).Raster = [];
        NeuronData{NeuronNo}.RasterPST(i).PST = [];
        NeuronData{NeuronNo}.RasterPST(i).Edges = Edges;

        TrialIndex = 0;

        % Find the number of trials with 'i' number of INs - both from the
        % beginning of bouts and within bouts

        Trials = length(find((Neural_INR{NeuronNo}.NoofINs == i))) + length(find(Neural_INR{NeuronNo}.WithinBoutNoofINs(:,1) == i));

        if (Trials > MinNumber)

            % Now first get all the data associated with the beginning of bouts

            BegTrials = find(Neural_INR{NeuronNo}.NoofINs == i);
            for j = 1:length(BegTrials),
                TrialSpikeTimes = [];
                TrialIndex = TrialIndex + 1;
                BoutSpikeTimes = Neural_INR{NeuronNo}.BoutDetails(BegTrials(j)).SpikeTimes + PreMotorLatency;
                INs = Neural_INR{NeuronNo}.INs{BegTrials(j)};
                % First get pre-time spikes
                PreTimeEnd = -(sum(MedianSyllDurs(1:i)) + sum(MedianGapDurs(1:i)));
                PreTimeStart = PreTimeEnd - PreTime;

                % First get the spikes occuring in the short pre-time period
                % before the first intro note
                TempSpikeIndices = find((BoutSpikeTimes >= (Neural_INR{NeuronNo}.BoutDetails(BegTrials(j)).onsets(INs(1)) - PreTime)) & (BoutSpikeTimes < Neural_INR{NeuronNo}.BoutDetails(BegTrials(j)).onsets(INs(1))));
                TempSpikeIndices = TempSpikeIndices(:);

                if (~isempty(TempSpikeIndices))
                    TrialSpikeTimes = [TrialSpikeTimes; (BoutSpikeTimes(TempSpikeIndices) - (Neural_INR{NeuronNo}.BoutDetails(BegTrials(j)).onsets(INs(1)) - PreTime) + PreTimeStart)];
                end

                TrialSpikeTimes = TrialSpikeTimes(:);
                % Now get the spikes during each IN and during each gap and
                % warp it to the median IN length and median gap length -
                % starting at the end, so that the end is the last IN

                for k = length(INs):-1:1,
                    % First get spikes associated with the gap
                    NextSyllONTime = Neural_INR{NeuronNo}.BoutDetails(BegTrials(j)).onsets(INs(k) + 1);
                    SyllOFFTime = Neural_INR{NeuronNo}.BoutDetails(BegTrials(j)).offsets(INs(k));
                    TempSpikeIndices = find((BoutSpikeTimes >= SyllOFFTime) & (BoutSpikeTimes < NextSyllONTime));
                    TempSpikeIndices = TempSpikeIndices(:);
                    if (~isempty(TempSpikeIndices))
                        WarpedSpikeTimes = (BoutSpikeTimes(TempSpikeIndices) - SyllOFFTime) * MedianGapDurs(length(INs) - k + 1) / (NextSyllONTime - SyllOFFTime);
                        WarpedSpikeTimes = WarpedSpikeTimes(:);
                        TrialSpikeTimes = [TrialSpikeTimes; (WarpedSpikeTimes - sum(MedianGapDurs(1:(length(INs) - k + 1))) - sum(MedianSyllDurs(1:(length(INs) - k))))];
                    end

                    % Now get the spikes associated with the IN
                    SyllONTime = Neural_INR{NeuronNo}.BoutDetails(BegTrials(j)).onsets(INs(k));
                    SyllOFFTime = Neural_INR{NeuronNo}.BoutDetails(BegTrials(j)).offsets(INs(k));
                    TempSpikeIndices = find((BoutSpikeTimes >= SyllONTime) & (BoutSpikeTimes < SyllOFFTime));
                    TempSpikeIndices = TempSpikeIndices(:);
                    if (~isempty(TempSpikeIndices))
                        WarpedSpikeTimes = (BoutSpikeTimes(TempSpikeIndices) - SyllONTime) * MedianSyllDurs(length(INs) - k + 1) / (SyllOFFTime - SyllONTime);
                        WarpedSpikeTimes = WarpedSpikeTimes(:);
                        TrialSpikeTimes = [TrialSpikeTimes; (WarpedSpikeTimes - sum(MedianGapDurs(1:(length(INs) - k + 1))) - sum(MedianSyllDurs(1:(length(INs) - k + 1))))];
                    end

                end
                % Now get the spikes associated with the first syllable
                SyllONTime = Neural_INR{NeuronNo}.BoutDetails(BegTrials(j)).onsets(INs(end) + 1);
                SyllOFFTime = Neural_INR{NeuronNo}.BoutDetails(BegTrials(j)).offsets(INs(end) + 1);
                TempSpikeIndices = find((BoutSpikeTimes >= SyllONTime) & (BoutSpikeTimes < SyllOFFTime));
                TempSpikeIndices = TempSpikeIndices(:);
                if (~isempty(TempSpikeIndices))
                    WarpedSpikeTimes = (BoutSpikeTimes(TempSpikeIndices) - SyllONTime) * MedianFirstSyllDur / (SyllOFFTime - SyllONTime);
                    WarpedSpikeTimes = WarpedSpikeTimes(:);
                    TrialSpikeTimes = [TrialSpikeTimes; WarpedSpikeTimes];
                end

                TrialSpikeTimes = sort(TrialSpikeTimes);

                NeuronData{NeuronNo}.RasterPST(i).TrialSpikeTimes{TrialIndex} = TrialSpikeTimes;
                NeuronData{NeuronNo}.RasterPST(i).Raster = [NeuronData{NeuronNo}.RasterPST(i).Raster; [TrialSpikeTimes ones(size(TrialSpikeTimes))*TrialIndex]];
                if (isempty(TrialSpikeTimes))
                    NeuronData{NeuronNo}.RasterPST(i).PST(TrialIndex, :) = zeros(size(Edges));
                else
                    NeuronData{NeuronNo}.RasterPST(i).PST(TrialIndex, :) = histc(TrialSpikeTimes, Edges)/BinSize;
                end

            end

            % Now get the data associated with within bouts

            WithinBoutTrials = find(Neural_INR{NeuronNo}.WithinBoutNoofINs(:,1) == i);
            for j = 1:length(WithinBoutTrials),
                TrialSpikeTimes = [];

                TrialIndex = TrialIndex + 1;
                BoutIndex = Neural_INR{NeuronNo}.WithinBoutINBoutIndices(WithinBoutTrials(j));
                BoutSpikeTimes = Neural_INR{NeuronNo}.BoutDetails(BoutIndex).SpikeTimes + PreMotorLatency;
                INs = Neural_INR{NeuronNo}.WithinBoutINs{WithinBoutTrials(j)};

                % First get pre-time spikes
                PreTimeEnd = -(sum(MedianSyllDurs(1:i)) + sum(MedianGapDurs(1:i)));
                PreTimeStart = PreTimeEnd - PreTime;

                % First get the spikes occuring in the short pre-time period
                % before the first intro note
                TempSpikeIndices = find((BoutSpikeTimes >= (Neural_INR{NeuronNo}.BoutDetails(BoutIndex).onsets(INs(1)) - PreTime)) & (BoutSpikeTimes < Neural_INR{NeuronNo}.BoutDetails(BoutIndex).onsets(INs(1))));
                TempSpikeIndices = TempSpikeIndices(:);
                if (~isempty(TempSpikeIndices))
                    TrialSpikeTimes = [TrialSpikeTimes; (BoutSpikeTimes(TempSpikeIndices) - (Neural_INR{NeuronNo}.BoutDetails(BoutIndex).onsets(INs(1)) - PreTime) + PreTimeStart)];
                end

                TrialSpikeTimes = TrialSpikeTimes(:);
                % Now get the spikes during each IN and during each gap and
                % warp it to the median IN length and median gap length -
                % starting at the end, so that the end is the last IN

                for k = length(INs):-1:1,
                    % First get spikes associated with the gap
                    NextSyllONTime = Neural_INR{NeuronNo}.BoutDetails(BoutIndex).onsets(INs(k) + 1);
                    SyllOFFTime = Neural_INR{NeuronNo}.BoutDetails(BoutIndex).offsets(INs(k));
                    TempSpikeIndices = find((BoutSpikeTimes >= SyllOFFTime) & (BoutSpikeTimes < NextSyllONTime));
                    TempSpikeIndices = TempSpikeIndices(:);
                    if (~isempty(TempSpikeIndices))
                        WarpedSpikeTimes = (BoutSpikeTimes(TempSpikeIndices) - SyllOFFTime) * MedianGapDurs(length(INs) - k + 1) / (NextSyllONTime - SyllOFFTime);
                        WarpedSpikeTimes = WarpedSpikeTimes(:);
                        TrialSpikeTimes = [TrialSpikeTimes; (WarpedSpikeTimes - sum(MedianGapDurs(1:(length(INs) - k + 1))) - sum(MedianSyllDurs(1:(length(INs) - k))))];
                    end

                    % Now get the spikes associated with the IN
                    SyllONTime = Neural_INR{NeuronNo}.BoutDetails(BoutIndex).onsets(INs(k));
                    SyllOFFTime = Neural_INR{NeuronNo}.BoutDetails(BoutIndex).offsets(INs(k));
                    TempSpikeIndices = find((BoutSpikeTimes >= SyllONTime) & (BoutSpikeTimes < SyllOFFTime));
                    TempSpikeIndices = TempSpikeIndices(:);
                    if (~isempty(TempSpikeIndices))
                        WarpedSpikeTimes = (BoutSpikeTimes(TempSpikeIndices) - SyllONTime) * MedianSyllDurs(length(INs) - k + 1) / (SyllOFFTime - SyllONTime);
                        WarpedSpikeTimes = WarpedSpikeTimes(:);
                        TrialSpikeTimes = [TrialSpikeTimes; (WarpedSpikeTimes - sum(MedianGapDurs(1:(length(INs) - k + 1))) - sum(MedianSyllDurs(1:(length(INs) - k + 1))))];
                    end
                end

                % Now get the spikes associated with the first syllable
                SyllONTime = Neural_INR{NeuronNo}.BoutDetails(BoutIndex).onsets(INs(end) + 1);
                SyllOFFTime = Neural_INR{NeuronNo}.BoutDetails(BoutIndex).offsets(INs(end) + 1);
                TempSpikeIndices = find((BoutSpikeTimes >= SyllONTime) & (BoutSpikeTimes < SyllOFFTime));
                TempSpikeIndices = TempSpikeIndices(:);
                if (~isempty(TempSpikeIndices))
                    WarpedSpikeTimes = (BoutSpikeTimes(TempSpikeIndices) - SyllONTime) * MedianFirstSyllDur / (SyllOFFTime - SyllONTime);
                    WarpedSpikeTimes = WarpedSpikeTimes(:);
                    TrialSpikeTimes = [TrialSpikeTimes; WarpedSpikeTimes];
                end

                TrialSpikeTimes = sort(TrialSpikeTimes);

                NeuronData{NeuronNo}.RasterPST(i).TrialSpikeTimes{TrialIndex} = TrialSpikeTimes;
                NeuronData{NeuronNo}.RasterPST(i).Raster = [NeuronData{NeuronNo}.RasterPST(i).Raster; [TrialSpikeTimes ones(size(TrialSpikeTimes))*TrialIndex]];
                if (isempty(TrialSpikeTimes))
                    NeuronData{NeuronNo}.RasterPST(i).PST(TrialIndex, :) = zeros(size(Edges));
                else
                    NeuronData{NeuronNo}.RasterPST(i).PST(TrialIndex, :) = histc(TrialSpikeTimes, Edges)/BinSize;
                end
            end
        end
    end
end

Colours = 'rbk';
figure(13);
RasterIncrement = 0;
ValidNums = [];
MeanPST = [];
for i = 1:length(NeuronData),
    if (size(NeuronData{i}.RasterPST(INNum).PST, 1) >= 2)
        ValidNums = [ValidNums; i];
        subplot(6,1,[1:4]);
        PlotRaster(NeuronData{i}.RasterPST(INNum).Raster, Colours(mod(i-1,length(Colours)) + 1), 0.25, RasterIncrement, min([size(NeuronData{i}.RasterPST(INNum).PST,1) MinPlotNumber]));
        hold on;
        RasterIncrement = RasterIncrement + min([size(NeuronData{i}.RasterPST(INNum).PST,1) MinPlotNumber]) + 1;
        
        subplot(6,1,[5:6]);
        hold on;
        % plot(NeuronData{i}.RasterPST(INNum).Edges, conv(mean(NeuronData{i}.RasterPST(INNum).PST), GaussWin, 'same'), Colours(mod(i-1, length(Colours)) + 1));
        MeanPST = [MeanPST; conv(mean(NeuronData{i}.RasterPST(INNum).PST), GaussWin, 'same')];
    end
end
subplot(6,1,[1:4]);
axis tight;
temp = axis;
axis([NeuronData{1}.RasterPST(INNum).Edges(1) NeuronData{1}.RasterPST(INNum).Edges(end) 0 1.05*temp(4)]);

for i = 1:INNum,
    YValues = [0.5 0.5 temp(4) temp(4)];
    fill([-(sum(MedianSyllDurs(1:i)) + sum(MedianGapDurs(1:i))) -(sum(MedianGapDurs(1:i)) + sum(MedianSyllDurs(1:i-1))) fliplr([-(sum(MedianSyllDurs(1:i)) + sum(MedianGapDurs(1:i))) -(sum(MedianGapDurs(1:i)) + sum(MedianSyllDurs(1:i-1)))])], YValues, 'k', 'EdgeColor', 'none', 'FaceColor', [0.75 0.75 0.75], 'FaceAlpha', 0.5);
end
YValues = [0.5 0.5 temp(4) temp(4)];
fill([0 MedianFirstSyllDur MedianFirstSyllDur 0], YValues, 'k', 'EdgeColor', 'none', 'FaceColor', [0.75 0.75 0.75], 'FaceAlpha', 0.5);

subplot(6,1,[5:6]);
hold on;
fill([NeuronData{1}.RasterPST(INNum).Edges fliplr(NeuronData{1}.RasterPST(INNum).Edges)], [(mean(MeanPST) + std(MeanPST)/sqrt(size(MeanPST,1))) fliplr((mean(MeanPST) - std(MeanPST)/sqrt(size(MeanPST,1))))], 'k', 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(NeuronData{1}.RasterPST(INNum).Edges, mean(MeanPST), 'k');
axis tight;
temp = axis;
axis([NeuronData{1}.RasterPST(INNum).Edges(1) NeuronData{1}.RasterPST(INNum).Edges(end) 0 1.05*temp(4)]);
for i = 1:INNum,
    fill([-(sum(MedianSyllDurs(1:i)) + sum(MedianGapDurs(1:i))) -(sum(MedianGapDurs(1:i)) + sum(MedianSyllDurs(1:i-1))) fliplr([-(sum(MedianSyllDurs(1:i)) + sum(MedianGapDurs(1:i))) -(sum(MedianGapDurs(1:i)) + sum(MedianSyllDurs(1:i-1)))])], [0.5 0.5 temp(4) temp(4)], 'k', 'EdgeColor', 'none', 'FaceColor', [0.75 0.75 0.75], 'FaceAlpha', 0.5);
end
YValues = [0.5 0.5 temp(4) temp(4)];
fill([0 MedianFirstSyllDur MedianFirstSyllDur 0], YValues, 'k', 'EdgeColor', 'none', 'FaceColor', [0.75 0.75 0.75], 'FaceAlpha', 0.5);

BinSize = 0.01;
Edges = [];
Edges = [Edges linspace(0, MedianFirstSyllDur, round(MedianFirstSyllDur/BinSize))];
for i = 1:INNum,
    StartPoint = -sum(MedianGapDurs(1:i)) - sum(MedianSyllDurs(1:i-1));
    EndPoint = -sum(MedianGapDurs(1:i-1)) - sum(MedianSyllDurs(1:i-1));
    Edges = [Edges linspace(StartPoint, EndPoint, round((EndPoint - StartPoint)/BinSize))];
    StartPoint = -sum(MedianGapDurs(1:i)) - sum(MedianSyllDurs(1:i));
    EndPoint = -sum(MedianGapDurs(1:i)) - sum(MedianSyllDurs(1:i-1));
    Edges = [Edges linspace(StartPoint, EndPoint, round((EndPoint - StartPoint)/BinSize))];
end
Edges = sort(Edges);
Edges = [Edges linspace((Edges(1) - 0.2), Edges(1), 0.2/BinSize)];
Edges = sort(Edges);

for i = 1:size(MeanPST,1),
    MeanSmoothPST(i,:) = interp1(NeuronData{1}.RasterPST(INNum).Edges, MeanPST(i,:), Edges);
end

if (~exist('PC', 'var'))
    [PC, Var] = eig(cov(MeanPST'));
end

Score = MeanPST'*PC;

if (~exist('Symbol', 'var'))
    Symbol = 'o';
end

figure(12);
hold on;
Start = 1;
End = sum(MedianGapDurs(1:INNum)) + sum(MedianSyllDurs(1:INNum));
End = find(NeuronData{1}.RasterPST(INNum).Edges > -End, 1, 'first');
plot3(Score(Start:End,1), Score(Start:End, 2), Score(Start:End,3), ['k', Symbol, '-']);

Colours = ['rgbcmk'];
Index = 1;
for i = INNum:-1:1,
    Index = Index + 1;
    figure(1);
    hold on;
    Start = End;
    End = sum(MedianGapDurs(1:i)) + sum(MedianSyllDurs(1:i-1));
    End = find(NeuronData{1}.RasterPST(INNum).Edges > -End, 1, 'first');
    plot3(Score(Start:End,1), Score(Start:End, 2), Score(Start:End,3), [Colours(mod(i-1, length(Colours))+1), 'o', '-']);
    
    Start = End;
    End = sum(MedianGapDurs(1:i-1)) + sum(MedianSyllDurs(1:i-1));
    End = find(NeuronData{1}.RasterPST(INNum).Edges > -End, 1, 'first');
    plot3(Score(Start:End,1), Score(Start:End, 2), Score(Start:End,3), [Colours(mod(i-1, length(Colours))+1), '+', '-']);
end

figure(11);
hold on;
Start = End;
plot3(Score(Start:end,1), Score(Start:end, 2), Score(Start:end,3), ['m', Symbol, '-']);

disp('Finished');