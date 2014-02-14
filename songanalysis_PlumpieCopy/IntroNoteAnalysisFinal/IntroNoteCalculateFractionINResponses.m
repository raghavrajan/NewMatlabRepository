function [FractionResponses] = IntroNoteCalculateFractionINResponses(Neural_INR, BaselineFR, PreMotorLatency, ThreshMultiplier)

FR = []; % has three columns - first one is firing rate during the intro note, second column is firing rate during the gap, third column considers gap and IN together
INIndex = 0;
for i = 1:length(Neural_INR.NoofINs),
    if (Neural_INR.NoofINs(i) > 0)
        BoutSpikeTimes = Neural_INR.BoutDetails(i).SpikeTimes;
        INs = Neural_INR.INs{i};
        for j = 1:length(INs),
            INIndex = INIndex + 1;
            Position(INIndex) = j - length(Neural_INR.NoofINs(i)) - 1;
            INOnset = Neural_INR.BoutDetails(i).onsets(INs(j)) - PreMotorLatency;
            INOffset = Neural_INR.BoutDetails(i).offsets(INs(j)) - PreMotorLatency;
            GapOnset = INOffset;
            GapOffset = Neural_INR.BoutDetails(i).onsets(INs(j)+1) - PreMotorLatency;
            
            FR(INIndex,1) = length(find((BoutSpikeTimes >= INOnset) & (BoutSpikeTimes < INOffset)))/(INOffset - INOnset);
            FR(INIndex,2) = length(find((BoutSpikeTimes >= GapOnset) & (BoutSpikeTimes < GapOffset)))/(GapOffset - GapOnset);
            FR(INIndex,3) = length(find((BoutSpikeTimes >= INOnset) & (BoutSpikeTimes < GapOffset)))/(GapOffset - INOnset);
        end
    end
end

for i = 1:size(Neural_INR.WithinBoutNoofINs, 1),
    if (Neural_INR.WithinBoutNoofINs(i,1) > 0)
        BoutSpikeTimes = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes;
        INs = Neural_INR.WithinBoutINs{i};
        for j = 1:length(INs),
            INIndex = INIndex + 1;
            Position(INIndex) = j - length(Neural_INR.NoofINs(Neural_INR.WithinBoutINBoutIndices(i))) - 1;
            INOnset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j)) - PreMotorLatency;
            INOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(INs(j)) - PreMotorLatency;
            GapOnset = INOffset;
            GapOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j)+1) - PreMotorLatency;
            
            FR(INIndex,1) = length(find((BoutSpikeTimes >= INOnset) & (BoutSpikeTimes < INOffset)))/(INOffset - INOnset);
            FR(INIndex,2) = length(find((BoutSpikeTimes >= GapOnset) & (BoutSpikeTimes < GapOffset)))/(GapOffset - GapOnset);
            FR(INIndex,3) = length(find((BoutSpikeTimes >= INOnset) & (BoutSpikeTimes < GapOffset)))/(GapOffset - INOnset);
        end
    end
end

UpperThreshold = mean(BaselineFR(:,1)) + ThreshMultiplier*std(BaselineFR(:,1));
LowerThreshold = mean(BaselineFR(:,1)) - ThreshMultiplier*std(BaselineFR(:,1));

for i = 1:size(FR,2),
    FractionResponses(i) = (length(find((FR(:,i) > UpperThreshold) | (FR(:,i) < LowerThreshold)))/length(FR(:,i)));
end