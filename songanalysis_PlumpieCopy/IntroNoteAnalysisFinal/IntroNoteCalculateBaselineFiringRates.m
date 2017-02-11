function [BaselineFR] = IntroNoteCalculateBaselineFiringRates(Neural_INR, BaselineDur)

for i = 1:length(Neural_INR.BoutDetails),
    BaselineOnsetSpikeIndices = find((Neural_INR.BoutDetails(i).SpikeTimes >= Neural_INR.BoutDetails(i).BoutOnset) & (Neural_INR.BoutDetails(i).SpikeTimes < (Neural_INR.BoutDetails(i).BoutOnset + BaselineDur)));
    BaselineFR(i,1) = length(BaselineOnsetSpikeIndices)/BaselineDur;
    
    IFRFs = 1/(Neural_INR.BoutDetails(i).IFR(1,2) - Neural_INR.BoutDetails(i).IFR(1,1));
    BaselineFR(i,3) = mean(Neural_INR.BoutDetails(i).IFR(2,1:round(BaselineDur*IFRFs)));
    
    BaselineOffsetSpikeIndices = find((Neural_INR.BoutDetails(i).SpikeTimes >= (Neural_INR.BoutDetails(i).BoutOffset - BaselineDur)) & (Neural_INR.BoutDetails(i).SpikeTimes < (Neural_INR.BoutDetails(i).BoutOffset)));
    BaselineFR(i,2) = length(BaselineOffsetSpikeIndices)/BaselineDur;
end
