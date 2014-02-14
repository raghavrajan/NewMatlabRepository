function [INFR, INEdges, MeanINFR, STDINFR, GapFR, GapEdges, MeanGapFR, STDGapFR, Position] = IntroNoteCalculateAverageINGapFRHist(Neural_INR, PreTime, BinSize)

GapPreMotorLag = 0.045;

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
            
            INPSTEdges = -PreTime:BinSize:0.2;
            Edges = -PreTime:BinSize:(GapOffset - INOnset - GapPreMotorLag);
            
            INFR(INIndex,:) = histc(BoutSpikeTimes, INPSTEdges + INOnset)/BinSize;
            INSpikeTimes{INIndex} = BoutSpikeTimes(find((BoutSpikeTimes >= (INPSTEdges(1) + INOnset)) & (BoutSpikeTimes < (INPSTEdges(end) + INOnset)))) - INOnset;
            
            GapEdges = -PreTime:BinSize:(GapOffset - GapOnset - GapPreMotorLag);
            GapSpikeTimes{INIndex} = BoutSpikeTimes(find((BoutSpikeTimes >= (GapEdges(1) + GapOnset)) & (BoutSpikeTimes < (GapEdges(end) + GapOnset)))) - GapOnset;
            GapDur(INIndex,1) = GapOffset - GapOnset + PreTime - GapPreMotorLag;
            
            GapFR{INIndex} = histc(BoutSpikeTimes, GapEdges + GapOnset)/BinSize;
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
                        
            INPSTEdges = -PreTime:BinSize:0.2;
            Edges = -PreTime:BinSize:(GapOffset - INOnset - GapPreMotorLag);
            INFR(INIndex,:) = histc(BoutSpikeTimes, INPSTEdges + INOnset)/BinSize;
            INSpikeTimes{INIndex} = BoutSpikeTimes(find((BoutSpikeTimes >= (INPSTEdges(1) + INOnset)) & (BoutSpikeTimes < (INPSTEdges(end) + INOnset)))) - INOnset;
            
            GapEdges = -PreTime:BinSize:(GapOffset - GapOnset - GapPreMotorLag);
            GapSpikeTimes{INIndex} = BoutSpikeTimes(find((BoutSpikeTimes >= (GapEdges(1) + GapOnset)) & (BoutSpikeTimes < (GapEdges(end) + GapOnset)))) - GapOnset;
            GapDur(INIndex,1) = GapOffset - GapOnset + PreTime - GapPreMotorLag;
            
            GapFR{INIndex} = histc(BoutSpikeTimes, GapEdges + GapOnset)/BinSize;
        end
    end
end
%INFR(find(Position(:,1) ~= -2),:) = [];

MeanINFR(1,:) = mean(INFR(:,1:end-1));
STDINFR(1,:) = std(INFR(:,1:end-1));

Edges = -PreTime:BinSize:10;
Edges = Edges + BinSize/2;
INEdges = INPSTEdges(1:end-1) + BinSize/2;

%GapFR(Position(:,1) ~= -2) = [];

Lens = cellfun(@length, GapFR);

ANOVAPST = [];
for i = 1:max(Lens),
    Indices = find(Lens >= i);
    if (length(Indices) >= 3)
        TempPST = [];
        for j = 1:length(Indices),
            TempPST = [TempPST; GapFR{Indices(j)}(i)];
        end
        ANOVAPST = [ANOVAPST; [TempPST ones(size(TempPST))*i]];
        MeanGapFR(1,i) = mean(TempPST);
        STDGapFR(1,i) = std(TempPST);
    end
end


Edges = -PreTime:BinSize:10;
Edges = Edges + BinSize/2;
GapEdges = Edges(1:length(MeanGapFR));

% for i = 1:1000,
%     clear TempINFR;
%     for j = 1:length(INSpikeTimes),
%         TempSpikeTimes = INSpikeTimes{j} + rand*(INPSTEdges(end) + PreTime);
%         TempSpikeTimes(find(TempSpikeTimes > (INPSTEdges(end)))) = TempSpikeTimes(find(TempSpikeTimes > (INPSTEdges(end)))) - (INPSTEdges(end) + PreTime);
%         TempINFR(j,:) = histc(TempSpikeTimes, INPSTEdges)/BinSize;
%     end
%     MeanTempINFR(1,:) = mean(TempINFR(:,1:end-1));
%     RandomINFRMaxMin(i,:) = [max(MeanTempINFR) min(MeanTempINFR)];
% end
% 
% RandomINFRMaxMin = sort(RandomINFRMaxMin);
% 
% figure;
% %subplot(1,2,1);
% plot(INEdges, MeanINFR, 'r');
% hold on;
% plot([INEdges(1) INEdges(end)], [RandomINFRMaxMin(round(0.95*size(RandomINFRMaxMin,1)),1) RandomINFRMaxMin(round(0.95*size(RandomINFRMaxMin,1)),1)], 'k--');
% plot([INEdges(1) INEdges(end)], [RandomINFRMaxMin(round(0.05*size(RandomINFRMaxMin,1)),2) RandomINFRMaxMin(round(0.05*size(RandomINFRMaxMin,1)),2)], 'k--');

% for k = 1:1000,
%     clear TempGapFR;
%     for j = 1:length(INSpikeTimes),
%         TempSpikeTimes = GapSpikeTimes{j} + rand*(GapDur(j));
%         TempSpikeTimes(find(TempSpikeTimes > (GapDur(j) - PreTime))) = TempSpikeTimes(find(TempSpikeTimes > (GapDur(j) - PreTime))) - GapDur(j);
%         TempGapEdges = -PreTime:BinSize:(GapDur(j) - PreTime);
%         TempGapFR{j} = histc(TempSpikeTimes, TempGapEdges)/BinSize;
%     end
%     Lens = cellfun(@length, TempGapFR);
% 
%     for i = 1:max(Lens),
%         Indices = find(Lens >= i);
%         if (length(Indices) >= 3)
%             TempPST = [];
%             for j = 1:length(Indices),
%                 TempPST = [TempPST; TempGapFR{Indices(j)}(i)];
%             end
%             MeanTempGapFR(1,i) = mean(TempPST);
%         end
%     end
%     RandomGapFRMaxMin(k,:) = [max(MeanTempGapFR) min(MeanTempGapFR)];
% end
% 
% RandomGapFRMaxMin = sort(RandomGapFRMaxMin);
% 
% subplot(1,2,2);
% plot(GapEdges, MeanGapFR);
% hold on;
% plot([GapEdges(1) GapEdges(end)], [RandomGapFRMaxMin(round(0.95*size(RandomGapFRMaxMin,1)),1) RandomGapFRMaxMin(round(0.95*size(RandomGapFRMaxMin,1)),1)], 'k--');
% plot([GapEdges(1) GapEdges(end)], [RandomGapFRMaxMin(round(0.05*size(RandomGapFRMaxMin,1)),2) RandomGapFRMaxMin(round(0.05*size(RandomGapFRMaxMin,1)),2)], 'k--');

disp('Finished');