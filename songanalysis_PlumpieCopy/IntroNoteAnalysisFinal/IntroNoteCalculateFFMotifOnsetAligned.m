function [FF, Edges] = IntroNoteCalculateFFMotifOnsetAligned(Neural_INR, PreTime, BinSize, AlignmentPoint)

GapPreMotorLag = 0.045;

Width = 0.005;
GaussianLen = 2;
IFRFs = 1/BinSize;
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (IFRFs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * IFRFs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * IFRFs) * (Width * IFRFs)));

Edges = -PreTime:0.001:0.2;

INIndex = 0;
for i = 1:length(Neural_INR.NoofINs),
    if (Neural_INR.NoofINs(i) > 0)
        INIndex = INIndex + 1;
        BoutSpikeTimes = Neural_INR.BoutDetails(i).SpikeTimes;
        if (AlignmentPoint ~= 1000)
            if (length(Neural_INR.INs{i}) < (abs(AlignmentPoint) + 1))
                INIndex = INIndex - 1;
                continue;
            else
                FFEdges = Edges + Neural_INR.BoutDetails(i).onsets(Neural_INR.INs{i}(end) + AlignmentPoint);
            end
        else
            FFEdges = Edges + Neural_INR.BoutDetails(i).onsets(Neural_INR.INs{i}(1));
            end 
        for j = 1:length(FFEdges),
            FF(INIndex,j) = length(find((BoutSpikeTimes >= (FFEdges(j))) & (BoutSpikeTimes < (FFEdges(j) + BinSize))));
        end
    end
end

for i = 1:size(Neural_INR.WithinBoutNoofINs, 1),
    if (Neural_INR.WithinBoutNoofINs(i,1) > 0)
        INIndex = INIndex + 1;
        BoutSpikeTimes = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes;
        if (AlignmentPoint ~= 1000)
            if (length(Neural_INR.WithinBoutINs{i}) < (abs(AlignmentPoint) + 1))
                INIndex = INIndex - 1;
                continue;
            else
                FFEdges = Edges + Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(Neural_INR.WithinBoutINs{i}(end) + AlignmentPoint);
            end
        else
            FFEdges = Edges + Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(Neural_INR.WithinBoutINs{i}(1));
        end
        for j = 1:length(FFEdges),
            FF(INIndex,j) = length(find((BoutSpikeTimes >= (FFEdges(j))) & (BoutSpikeTimes < (FFEdges(j) + BinSize))));
        end
    end
end

FF = var(FF) - mean(FF);
disp('Finished');