function [INFR, GapFR, Edges, Position, IFREdges, INIFR, GapIFR, MeanGapIFR, MeanPosGapIFR] = IntroNoteCalculateAverageINGapFR(Neural_INR, PreTime, PostTime, BinSize)

GapPreMotorLag = 0.05;

Edges = -PreTime:BinSize:PostTime;
IFREdges = -PreTime:1/2000:PostTime;
IFRWiderEdges = [(IFREdges(1) - 1/2000) IFREdges (IFREdges(end) + 1/2000)];

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
            Position(INIndex,1) = j - length(INs) - 1;
            Position(INIndex,2) = j;
            INOnset = Neural_INR.BoutDetails(i).onsets(INs(j));
            INOffset = Neural_INR.BoutDetails(i).offsets(INs(j));
            GapOnset = INOffset;
            GapOffset = Neural_INR.BoutDetails(i).onsets(INs(j)+1);
            
            INFR(INIndex,:) = histc(BoutSpikeTimes(find((BoutSpikeTimes >= (INOnset + Edges(1))) & (BoutSpikeTimes < (INOnset + Edges(end))))), Edges + INOnset);
            GapFR(INIndex,:) = histc(BoutSpikeTimes(find((BoutSpikeTimes >= (GapOnset + Edges(1))) & (BoutSpikeTimes < (GapOnset + Edges(end))))), Edges + GapOnset);
            
            INIFRIndices = find((Neural_INR.BoutDetails(i).IFR(1,:) >= (INOnset + Edges(1))) & (Neural_INR.BoutDetails(i).IFR(1,:) <= (INOnset + Edges(end))));
            TempIFR = spline(Neural_INR.BoutDetails(i).IFR(1, INIFRIndices), Neural_INR.BoutDetails(i).IFR(2, INIFRIndices), INOnset + IFRWiderEdges);
            INIFR(INIndex,:) = TempIFR(2:end-1);
            
            GapIFRIndices = find((Neural_INR.BoutDetails(i).IFR(1,:) >= (INOnset + Edges(1))) & (Neural_INR.BoutDetails(i).IFR(1,:) <= (GapOffset - GapPreMotorLag)));
            GapEdges = Edges(1):1/2000:(GapOffset - INOnset - GapPreMotorLag);
            GapWiderEdges = [(GapEdges(1)-1/2000) GapEdges (GapEdges(end)+1/2000)];
            TempIFR = spline(Neural_INR.BoutDetails(i).IFR(1, GapIFRIndices), Neural_INR.BoutDetails(i).IFR(2, GapIFRIndices), GapWiderEdges + INOnset);            
            TempIFR = conv(TempIFR, GaussWin, 'same');
            GapIFR{INIndex} = [GapEdges; TempIFR(2:end-1)];
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
            INOnset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j));
            INOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(INs(j));
            GapOnset = INOffset;
            GapOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j)+1);
                        
            INFR(INIndex,:) = histc(BoutSpikeTimes(find((BoutSpikeTimes >= (INOnset + Edges(1))) & (BoutSpikeTimes < (INOnset + Edges(end))))), Edges + INOnset);
            GapFR(INIndex,:) = histc(BoutSpikeTimes(find((BoutSpikeTimes >= (GapOnset + Edges(1))) & (BoutSpikeTimes < (GapOnset + Edges(end))))), Edges + GapOnset);
            
            INIFRIndices = find((Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).IFR(1,:) >= (INOnset + Edges(1))) & (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).IFR(1,:) <= (INOnset + Edges(end))));
            TempIFR = spline(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).IFR(1, INIFRIndices), Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).IFR(2, INIFRIndices), INOnset + IFRWiderEdges);
            INIFR(INIndex,:) = TempIFR(2:end-1);
            
            GapIFRIndices = find((Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).IFR(1,:) >= (INOnset + Edges(1))) & (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).IFR(1,:) <= (GapOffset - GapPreMotorLag)));
            GapEdges = Edges(1):1/2000:(GapOffset - INOnset - GapPreMotorLag);
            GapWiderEdges = [(GapEdges(1)-1/2000) GapEdges (GapEdges(end)+1/2000)];
            TempIFR = spline(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).IFR(1, GapIFRIndices), Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).IFR(2, GapIFRIndices), GapWiderEdges + INOnset);            
            TempIFR = conv(TempIFR, GaussWin, 'same');
            GapIFR{INIndex} = [GapEdges; TempIFR(2:end-1)];            
        end
    end
end

INFR = INFR./BinSize;
GapFR = GapFR./BinSize;

INFR = INFR(:,1:end-1);
GapFR = GapFR(:,1:end-1);

Edges = Edges(1:end-1) + BinSize/2;

Lens = cellfun(@length, GapIFR);

for i = 1:max(Lens),
    Indices = find(Lens >= i);
    if (length(Indices) >= 3)
        TempPST = [];
        for j = 1:length(Indices),
        TempPST = [TempPST; GapIFR{Indices(j)}(2,i)];
        end
        TempMeanGapIFR(1,i) = mean(TempPST);
    end
end

[i, j] = max(Lens);
MeanGapIFR = [GapIFR{j}(1,1:length(TempMeanGapIFR)); TempMeanGapIFR];


% TempMeanGapIFR = [];
% LastGaps = find(Position(:,1) == -1);
% 
% Lens = cellfun(@length, GapIFR);
% for i = 1:max(Lens),
%     Indices = find(Lens(LastGaps) >= i);
%     if (length(Indices) >= 3)
%         Indices = find(Lens >= i);
%         TempPST = [];
%         for j = 1:length(Indices),
%             if (~isempty(find(LastGaps == Indices(j))))
%                 TempPST = [TempPST; GapIFR{Indices(j)}(2,i)];
%             end
%         end
%         TempMeanGapIFR(1,i) = mean(TempPST);
%     end
% end
% 
% clear MeanGapIFR;
% [i, j] = max(Lens);
% MeanGapIFR = [GapIFR{j}(1,1:length(TempMeanGapIFR)); TempMeanGapIFR];

disp('Finished');