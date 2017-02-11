function [] = IntroNotePopulationNeuralAnalysis(Neural_INR)

Motif = 'abcdecd';
PreMotorLag = 0.045;

SyllDur = [];
GapDur = [];
Index = 0;
for i = 1:length(Neural_INR),
    Trials = find(Neural_INR{i}.NoofINs == 4);
    for j = 1:length(Trials),
        if (length(Neural_INR{i}.BoutDetails(Trials(j)).labels) >= (Neural_INR{i}.INs{Trials(j)}(end) + 7))
            if (strfind(Neural_INR{i}.BoutDetails(Trials(j)).labels((Neural_INR{i}.INs{Trials(j)}(end) + 1):(Neural_INR{i}.INs{Trials(j)}(end) + 7)), Motif))
                Index = Index + 1;
                for k = Neural_INR{i}.INs{Trials(j)}(1):(Neural_INR{i}.INs{Trials(j)}(end) + 7),
                    SyllDur(Index,(k - Neural_INR{i}.INs{Trials(j)}(1) + 1)) = Neural_INR{i}.BoutDetails(Trials(j)).offsets(k) - Neural_INR{i}.BoutDetails(Trials(j)).onsets(k);
                    if (k ~= (Neural_INR{i}.INs{Trials(j)}(end) + 7))
                        GapDur(Index,(k - Neural_INR{i}.INs{Trials(j)}(1) + 1)) = Neural_INR{i}.BoutDetails(Trials(j)).onsets(k + 1) - Neural_INR{i}.BoutDetails(Trials(j)).offsets(k);
                    end
                end
            end
        end
    end
end

MedianSyllDurs = median(SyllDur);
MedianGapDurs = median(GapDur);

Edges = 0:0.005:(sum(MedianGapDurs) + sum(MedianSyllDurs));
Width = 0.005;
GaussianLen = 4;
IFRFs = 1/0.005;
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (IFRFs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * IFRFs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * IFRFs) * (Width * IFRFs)));

for i = 1:length(Neural_INR),
    Trials = find(Neural_INR{i}.NoofINs == 4);
    Index = 0;
    for j = 1:length(Trials),
        if (length(Neural_INR{i}.BoutDetails(Trials(j)).labels) >= (Neural_INR{i}.INs{Trials(j)}(end) + 7))
            if (strfind(Neural_INR{i}.BoutDetails(Trials(j)).labels((Neural_INR{i}.INs{Trials(j)}(end) + 1):(Neural_INR{i}.INs{Trials(j)}(end) + 7)), Motif))
                Index = Index + 1;
                TempSpikeTimes = [];
                BoutSpikeTimes = Neural_INR{i}.BoutDetails(Trials(j)).SpikeTimes;
                for k = Neural_INR{i}.INs{Trials(j)}(1):(Neural_INR{i}.INs{Trials(j)}(end) + 7),
                    SyllIndex = k - Neural_INR{i}.INs{Trials(j)}(1) + 1;
                    SyllDur = Neural_INR{i}.BoutDetails(Trials(j)).offsets(k) - Neural_INR{i}.BoutDetails(Trials(j)).onsets(k);
                    SyllSpikeIndices = find((BoutSpikeTimes >= (Neural_INR{i}.BoutDetails(Trials(j)).onsets(k) - PreMotorLag)) & (BoutSpikeTimes < (Neural_INR{i}.BoutDetails(Trials(j)).offsets(k) - PreMotorLag)));
                    if (~isempty(SyllSpikeIndices))
                        SyllSpikes = (BoutSpikeTimes(SyllSpikeIndices) - (Neural_INR{i}.BoutDetails(Trials(j)).onsets(k) - PreMotorLag)) * MedianSyllDurs(SyllIndex)/SyllDur;
                        SyllSpikes = SyllSpikes(:) + sum(MedianGapDurs(1:SyllIndex-1)) + sum(MedianSyllDurs(1:SyllIndex-1));
                        TempSpikeTimes = [TempSpikeTimes; SyllSpikes];
                    end
                    if (k ~= (Neural_INR{i}.INs{Trials(j)}(end) + 7))
                        GapDur = Neural_INR{i}.BoutDetails(Trials(j)).onsets(k + 1) - Neural_INR{i}.BoutDetails(Trials(j)).offsets(k);
                        GapSpikeIndices = find((BoutSpikeTimes >= (Neural_INR{i}.BoutDetails(Trials(j)).offsets(k) - PreMotorLag)) & (BoutSpikeTimes < (Neural_INR{i}.BoutDetails(Trials(j)).onsets(k + 1) - PreMotorLag)));
                        if (~isempty(GapSpikeIndices))
                            GapSpikes = (BoutSpikeTimes(GapSpikeIndices) - (Neural_INR{i}.BoutDetails(Trials(j)).offsets(k) - PreMotorLag)) * MedianGapDurs(SyllIndex)/GapDur;
                            GapSpikes = GapSpikes(:) + sum(MedianGapDurs(1:SyllIndex-1)) + sum(MedianSyllDurs(1:SyllIndex));
                            TempSpikeTimes = [TempSpikeTimes; GapSpikes];
                        end
                    end
                end
                NeuronPST{i}(Index,:) = histc(TempSpikeTimes, Edges);
                SmoothNeuronPST{i}(Index,:) = conv(histc(TempSpikeTimes, Edges), GaussWin, 'same');
            end
        end
    end
end

disp('Finished');