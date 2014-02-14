function [INDetails] = IntroNoteAnalysisGetSyllGapSpikeTimes(Neural_INR, PreMotorLag, BinSize)

Index = 0;
IFRFs = 2000;
FirstSyllIndex = 0;

PreTime = 0.1;
PostTime = 0.1;
FullPSTEdges = -PreTime:BinSize:PostTime;

for i = 1:length(Neural_INR.NoofINs),
    if (Neural_INR.NoofINs(i) < 1)
        continue;
    end
    BoutSpikeTimes = Neural_INR.BoutDetails(i).SpikeTimes;
    BoutIFR = Neural_INR.BoutDetails(i).IFR;
    for j = 1:length(Neural_INR.INs{i}),
        Index = Index + 1;
        IN(Index).BoutBeginnings = 1;
        IN(Index).SequenceNo = i;
        IN(Index).Position = j - Neural_INR.NoofINs(i) - 1;
        IN(Index).PositionFromFirst = j;
        IN(Index).PreMotorLag = PreMotorLag;
        IN(Index).NoofINs = Neural_INR.NoofINs(i);
        
        BoutOnset = Neural_INR.BoutDetails(i).onsets(1) - 1.5;
        
        IN(Index).BaselineRate = length(find((BoutSpikeTimes >= BoutOnset) & (BoutSpikeTimes <= (BoutOnset + 0.3))))/0.3;
        
        IN(Index).INStartTime = Neural_INR.BoutDetails(i).onsets(Neural_INR.INs{i}(j));
        IN(Index).INEndTime = Neural_INR.BoutDetails(i).offsets(Neural_INR.INs{i}(j));
        IN(Index).INSpikeTimes = BoutSpikeTimes(find((BoutSpikeTimes >= (IN(Index).INStartTime - PreMotorLag)) & (BoutSpikeTimes <= (IN(Index).INEndTime - PreMotorLag)))) - (IN(Index).INStartTime - PreMotorLag);
        IN(Index).INFR = length(IN(Index).INSpikeTimes)/(IN(Index).INEndTime - IN(Index).INStartTime);
        
        IN(Index).INFullSpikeTimes = BoutSpikeTimes(find((BoutSpikeTimes >= (IN(Index).INStartTime - PreTime)) & (BoutSpikeTimes <= (IN(Index).INStartTime + PostTime)))) - (IN(Index).INStartTime);
        IN(Index).INFullPST = histc(IN(Index).INFullSpikeTimes, FullPSTEdges);
        IN(Index).INFullPST = IN(Index).INFullPST(:);
        
        IN(Index).INDur = IN(Index).INEndTime - IN(Index).INStartTime;
        IN(Index).INAmp = Neural_INR.BoutDetails(i).Feats(Neural_INR.INs{i}(j), 2);
        IN(Index).INEnt = Neural_INR.BoutDetails(i).Feats(Neural_INR.INs{i}(j), 3);
        IN(Index).INFreq = Neural_INR.BoutDetails(i).Feats(Neural_INR.INs{i}(j), 4);
        
        TempIFRIndices =  find((BoutIFR(1,:) >= (IN(Index).INStartTime - PreMotorLag)) & (BoutIFR(1,:) <= (IN(Index).INEndTime - PreMotorLag)));
        IN(Index).INIFR(1,:) = (IN(Index).INStartTime:1/IFRFs:IN(Index).INEndTime) - PreMotorLag;
        IN(Index).INIFR(2,:) = interp1(BoutIFR(1, TempIFRIndices), BoutIFR(2, TempIFRIndices), IN(Index).INIFR(1,:));
        IN(Index).INIFR(1,:) = IN(Index).INIFR(1,:) - (IN(Index).INStartTime - PreMotorLag);
        
        TempIFRIndices =  find((BoutIFR(1,:) >= (IN(Index).INStartTime - PreTime)) & (BoutIFR(1,:) <= (IN(Index).INEndTime - PostTime)));
        IN(Index).INFullIFR(1,:) = ((IN(Index).INStartTime - PreTime):1/IFRFs:(IN(Index).INStartTime + PostTime));
        IN(Index).INFullIFR(2,:) = interp1(BoutIFR(1, TempIFRIndices), BoutIFR(2, TempIFRIndices), IN(Index).INFullIFR(1,:));
        IN(Index).INFullIFR(1,:) = IN(Index).INFullIFR(1,:) - IN(Index).INStartTime;
        
        IN(Index).GapStartTime = Neural_INR.BoutDetails(i).offsets(Neural_INR.INs{i}(j));
        IN(Index).GapEndTime = Neural_INR.BoutDetails(i).onsets(Neural_INR.INs{i}(j) + 1);
        IN(Index).GapSpikeTimes = BoutSpikeTimes(find((BoutSpikeTimes >= (IN(Index).GapStartTime - PreMotorLag)) & (BoutSpikeTimes <= (IN(Index).GapEndTime - PreMotorLag)))) - (IN(Index).GapStartTime - PreMotorLag);
        IN(Index).GapDur = IN(Index).GapEndTime - (IN(Index).GapStartTime - PreMotorLag);
        IN(Index).GapFR = length(IN(Index).GapSpikeTimes)/(IN(Index).GapEndTime - IN(Index).GapStartTime);
        
        TempIFRIndices =  find((BoutIFR(1,:) >= (IN(Index).GapStartTime - PreMotorLag)) & (BoutIFR(1,:) <= (IN(Index).GapEndTime - PreMotorLag)));
        IN(Index).GapIFR(1,:) = (IN(Index).GapStartTime:1/IFRFs:IN(Index).GapEndTime) - PreMotorLag;
        IN(Index).GapIFR(2,:) = interp1(BoutIFR(1, TempIFRIndices), BoutIFR(2, TempIFRIndices), IN(Index).GapIFR(1,:));
        IN(Index).GapIFR(1,:) = IN(Index).GapIFR(1,:) - IN(Index).GapStartTime;
        
        IN(Index).INGapStartTime = Neural_INR.BoutDetails(i).onsets(Neural_INR.INs{i}(j));
        IN(Index).INGapEndTime = Neural_INR.BoutDetails(i).onsets(Neural_INR.INs{i}(j) + 1);
        IN(Index).INGapSpikeTimes = BoutSpikeTimes(find((BoutSpikeTimes >= (IN(Index).INGapStartTime - PreMotorLag)) & (BoutSpikeTimes <= (IN(Index).INGapEndTime - PreMotorLag)))) - (IN(Index).INGapStartTime - PreMotorLag);
        IN(Index).INGapDur = IN(Index).INGapEndTime - IN(Index).INGapStartTime;
        IN(Index).INGapFR = length(IN(Index).INGapSpikeTimes)/(IN(Index).INGapEndTime - IN(Index).INGapStartTime);
        
        TempIFRIndices =  find((BoutIFR(1,:) >= (IN(Index).INGapStartTime - PreMotorLag)) & (BoutIFR(1,:) <= (IN(Index).INGapEndTime - PreMotorLag)));
        IN(Index).INGapIFR(1,:) = (IN(Index).INGapStartTime:1/IFRFs:IN(Index).INGapEndTime) - PreMotorLag;
        IN(Index).INGapIFR(2,:) = interp1(BoutIFR(1, TempIFRIndices), BoutIFR(2, TempIFRIndices), IN(Index).INGapIFR(1,:));
        IN(Index).INGapIFR(1,:) = IN(Index).INGapIFR(1,:) - (IN(Index).INGapStartTime - PreMotorLag);
    end
    
    FirstSyllIndex = FirstSyllIndex + 1;
    FirstSyll(FirstSyllIndex).BoutBeginnings = 1;
    FirstSyll(FirstSyllIndex).SequenceNo = i;
    FirstSyll(FirstSyllIndex).Position = 0;
    FirstSyll(FirstSyllIndex).PreMotorLag = PreMotorLag;

    FirstSyll(FirstSyllIndex).FirstSyllStartTime = Neural_INR.BoutDetails(i).onsets(Neural_INR.INs{i}(end) + 1);
    FirstSyll(FirstSyllIndex).FirstSyllEndTime = Neural_INR.BoutDetails(i).offsets(Neural_INR.INs{i}(end) + 1);
    FirstSyll(FirstSyllIndex).FirstSyllSpikeTimes = BoutSpikeTimes(find((BoutSpikeTimes >= (FirstSyll(FirstSyllIndex).FirstSyllStartTime - PreMotorLag)) & (BoutSpikeTimes <= (FirstSyll(FirstSyllIndex).FirstSyllEndTime - PreMotorLag)))) - (FirstSyll(FirstSyllIndex).FirstSyllStartTime - PreMotorLag);
    FirstSyll(FirstSyllIndex).FirstSyllDur = FirstSyll(FirstSyllIndex).FirstSyllEndTime - FirstSyll(FirstSyllIndex).FirstSyllStartTime;
    FirstSyll(FirstSyllIndex).FirstSyllAmp = Neural_INR.BoutDetails(i).Feats(Neural_INR.INs{i}(end)+1, 2);
    FirstSyll(FirstSyllIndex).FirstSyllEnt = Neural_INR.BoutDetails(i).Feats(Neural_INR.INs{i}(end)+1, 3);
    FirstSyll(FirstSyllIndex).FirstSyllFreq = Neural_INR.BoutDetails(i).Feats(Neural_INR.INs{i}(end)+1, 4);
    
%     if (length(Neural_INR.BoutDetails(i).offsets) >= (Neural_INR.INs{i}(end) + 2))
%         FirstSyll(FirstSyllIndex).FirstSyllEndTime = Neural_INR.BoutDetails(i).offsets(Neural_INR.INs{i}(end) + 2);
%     else
%         FirstSyll(FirstSyllIndex).FirstSyllEndTime = Neural_INR.BoutDetails(i).offsets(Neural_INR.INs{i}(end) + 1);
%     end
%     FirstSyll(FirstSyllIndex).FirstSyllSpikeTimes = BoutSpikeTimes(find((BoutSpikeTimes >= (FirstSyll(FirstSyllIndex).FirstSyllStartTime - PreMotorLag)) & (BoutSpikeTimes <= (FirstSyll(FirstSyllIndex).FirstSyllEndTime - PreMotorLag)))) - (FirstSyll(FirstSyllIndex).FirstSyllStartTime - PreMotorLag);
%     FirstSyll(FirstSyllIndex).FirstSyllDur = FirstSyll(FirstSyllIndex).FirstSyllEndTime - FirstSyll(FirstSyllIndex).FirstSyllStartTime;
%     FirstSyll(FirstSyllIndex).FirstSyllAmp = Neural_INR.BoutDetails(i).Feats(Neural_INR.INs{i}(end)+1, 2);
%     FirstSyll(FirstSyllIndex).FirstSyllEnt = Neural_INR.BoutDetails(i).Feats(Neural_INR.INs{i}(end)+1, 3);
%     FirstSyll(FirstSyllIndex).FirstSyllFreq = Neural_INR.BoutDetails(i).Feats(Neural_INR.INs{i}(end)+1, 4);

    TempIFRIndices =  find((BoutIFR(1,:) >= (FirstSyll(FirstSyllIndex).FirstSyllStartTime - PreMotorLag)) & (BoutIFR(1,:) <= (FirstSyll(FirstSyllIndex).FirstSyllEndTime - PreMotorLag)));
    FirstSyll(FirstSyllIndex).FirstSyllIFR(1,:) = (FirstSyll(FirstSyllIndex).FirstSyllStartTime:1/IFRFs:FirstSyll(FirstSyllIndex).FirstSyllEndTime) - PreMotorLag;
    FirstSyll(FirstSyllIndex).FirstSyllIFR(2,:) = interp1(BoutIFR(1, TempIFRIndices), BoutIFR(2, TempIFRIndices), FirstSyll(FirstSyllIndex).FirstSyllIFR(1,:));
    FirstSyll(FirstSyllIndex).FirstSyllIFR(1,:) = FirstSyll(FirstSyllIndex).FirstSyllIFR(1,:) - (FirstSyll(FirstSyllIndex).FirstSyllStartTime - PreMotorLag);
end

for i = 1:size(Neural_INR.WithinBoutNoofINs, 1),
    if (Neural_INR.WithinBoutNoofINs(i,1) < 1)
        continue;
    end
    BoutNo = Neural_INR.WithinBoutINBoutIndices(i);
    BoutSpikeTimes = Neural_INR.BoutDetails(BoutNo).SpikeTimes;
    BoutIFR = Neural_INR.BoutDetails(BoutNo).IFR;
    for j = 1:length(Neural_INR.WithinBoutINs{i}),
        Index = Index + 1;
        IN(Index).BoutBeginnings = 0;
        IN(Index).SequenceNo = BoutNo;
        IN(Index).Position = j - Neural_INR.WithinBoutNoofINs(i, 1) - 1;
        IN(Index).PositionFromFirst = j;
        IN(Index).PreMotorLag = PreMotorLag;
        IN(Index).NoofINs = Neural_INR.WithinBoutNoofINs(i,1);
    
        IN(Index).INStartTime = Neural_INR.BoutDetails(BoutNo).onsets(Neural_INR.WithinBoutINs{i}(j));
        IN(Index).INEndTime = Neural_INR.BoutDetails(BoutNo).offsets(Neural_INR.WithinBoutINs{i}(j));
        IN(Index).INSpikeTimes = BoutSpikeTimes(find((BoutSpikeTimes >= (IN(Index).INStartTime - PreMotorLag)) & (BoutSpikeTimes <= (IN(Index).INEndTime - PreMotorLag)))) - (IN(Index).INStartTime - PreMotorLag);
        IN(Index).INFR = length(IN(Index).INSpikeTimes)/(IN(Index).INEndTime - IN(Index).INStartTime);
        
        IN(Index).INFullSpikeTimes = BoutSpikeTimes(find((BoutSpikeTimes >= (IN(Index).INStartTime - PreTime)) & (BoutSpikeTimes <= (IN(Index).INStartTime + PostTime)))) - (IN(Index).INStartTime);
        IN(Index).INFullPST = histc(IN(Index).INFullSpikeTimes, FullPSTEdges);
        IN(Index).INFullPST = IN(Index).INFullPST(:);

        IN(Index).INDur = IN(Index).INEndTime - IN(Index).INStartTime;
        IN(Index).INAmp = Neural_INR.BoutDetails(BoutNo).Feats(Neural_INR.WithinBoutINs{i}(j), 2);
        IN(Index).INEnt = Neural_INR.BoutDetails(BoutNo).Feats(Neural_INR.WithinBoutINs{i}(j), 3);
        IN(Index).INFreq = Neural_INR.BoutDetails(BoutNo).Feats(Neural_INR.WithinBoutINs{i}(j), 4);
        
        TempIFRIndices =  find((BoutIFR(1,:) >= (IN(Index).INStartTime - PreMotorLag)) & (BoutIFR(1,:) <= (IN(Index).INEndTime - PreMotorLag)));
        IN(Index).INIFR(1,:) = (IN(Index).INStartTime:1/IFRFs:IN(Index).INEndTime) - PreMotorLag;
        IN(Index).INIFR(2,:) = interp1(BoutIFR(1, TempIFRIndices), BoutIFR(2, TempIFRIndices), IN(Index).INIFR(1,:));
        IN(Index).INIFR(1,:) = IN(Index).INIFR(1,:) - (IN(Index).INStartTime - PreMotorLag);
        
        TempIFRIndices =  find((BoutIFR(1,:) >= (IN(Index).INStartTime - PreTime)) & (BoutIFR(1,:) <= (IN(Index).INStartTime + PostTime)));
        IN(Index).INFullIFR(1,:) = (IN(Index).INStartTime - PreTime):1/IFRFs:(IN(Index).INStartTime + PostTime);
        IN(Index).INFullIFR(2,:) = interp1(BoutIFR(1, TempIFRIndices), BoutIFR(2, TempIFRIndices), IN(Index).INFullIFR(1,:));
        IN(Index).INFullIFR(1,:) = IN(Index).INFullIFR(1,:) - IN(Index).INStartTime;
        
        IN(Index).GapStartTime = Neural_INR.BoutDetails(BoutNo).offsets(Neural_INR.WithinBoutINs{i}(j));
        IN(Index).GapEndTime = Neural_INR.BoutDetails(BoutNo).onsets(Neural_INR.WithinBoutINs{i}(j) + 1);
        IN(Index).GapSpikeTimes = BoutSpikeTimes(find((BoutSpikeTimes >= (IN(Index).GapStartTime - PreMotorLag)) & (BoutSpikeTimes <= (IN(Index).GapEndTime - PreMotorLag)))) - (IN(Index).GapStartTime - PreMotorLag);
        IN(Index).GapDur = IN(Index).GapEndTime - (IN(Index).GapStartTime - PreMotorLag);
        IN(Index).GapFR = length(IN(Index).GapSpikeTimes)/(IN(Index).GapEndTime - IN(Index).GapStartTime);
        
        TempIFRIndices =  find((BoutIFR(1,:) >= (IN(Index).GapStartTime - PreMotorLag)) & (BoutIFR(1,:) <= (IN(Index).GapEndTime - PreMotorLag)));
        IN(Index).GapIFR(1,:) = (IN(Index).GapStartTime:1/IFRFs:IN(Index).GapEndTime) - PreMotorLag;
        IN(Index).GapIFR(2,:) = interp1(BoutIFR(1, TempIFRIndices), BoutIFR(2, TempIFRIndices), IN(Index).GapIFR(1,:));
        IN(Index).GapIFR(1,:) = IN(Index).GapIFR(1,:) - IN(Index).GapStartTime;
        
        IN(Index).INGapStartTime = Neural_INR.BoutDetails(BoutNo).onsets(Neural_INR.WithinBoutINs{i}(j));
        IN(Index).INGapEndTime = Neural_INR.BoutDetails(BoutNo).onsets(Neural_INR.WithinBoutINs{i}(j) + 1);
        IN(Index).INGapSpikeTimes = BoutSpikeTimes(find((BoutSpikeTimes >= (IN(Index).INGapStartTime - PreMotorLag)) & (BoutSpikeTimes <= (IN(Index).INGapEndTime - PreMotorLag)))) - (IN(Index).INGapStartTime - PreMotorLag);
        IN(Index).INGapDur = IN(Index).INGapEndTime - IN(Index).INGapStartTime;
        IN(Index).INGapFR = length(IN(Index).INGapSpikeTimes)/(IN(Index).INGapEndTime - IN(Index).INGapStartTime);
        
        TempIFRIndices =  find((BoutIFR(1,:) >= (IN(Index).INGapStartTime - PreMotorLag)) & (BoutIFR(1,:) <= (IN(Index).INGapEndTime - PreMotorLag)));
        IN(Index).INGapIFR(1,:) = (IN(Index).INGapStartTime:1/IFRFs:IN(Index).INGapEndTime) - PreMotorLag;
        IN(Index).INGapIFR(2,:) = interp1(BoutIFR(1, TempIFRIndices), BoutIFR(2, TempIFRIndices), IN(Index).INGapIFR(1,:));
        IN(Index).INGapIFR(1,:) = IN(Index).INGapIFR(1,:) - (IN(Index).INGapStartTime - PreMotorLag);
    end
    
    FirstSyllIndex = FirstSyllIndex + 1;
    FirstSyll(FirstSyllIndex).BoutBeginnings = 0;
    FirstSyll(FirstSyllIndex).SequenceNo = BoutNo;
    FirstSyll(FirstSyllIndex).Position = 0;
    FirstSyll(FirstSyllIndex).PreMotorLag = PreMotorLag;

    FirstSyll(FirstSyllIndex).FirstSyllStartTime = Neural_INR.BoutDetails(BoutNo).onsets(Neural_INR.WithinBoutINs{i}(end) + 1);
    if (length(Neural_INR.BoutDetails(BoutNo).offsets) >= (Neural_INR.WithinBoutINs{i}(end) + 2))
        FirstSyll(FirstSyllIndex).FirstSyllEndTime = Neural_INR.BoutDetails(BoutNo).offsets(Neural_INR.WithinBoutINs{i}(end) + 2);
    else
        FirstSyll(FirstSyllIndex).FirstSyllEndTime = Neural_INR.BoutDetails(BoutNo).offsets(Neural_INR.WithinBoutINs{i}(end) + 1);
    end
    FirstSyll(FirstSyllIndex).FirstSyllSpikeTimes = BoutSpikeTimes(find((BoutSpikeTimes >= (FirstSyll(FirstSyllIndex).FirstSyllStartTime - PreMotorLag)) & (BoutSpikeTimes <= (FirstSyll(FirstSyllIndex).FirstSyllEndTime - PreMotorLag)))) - (FirstSyll(FirstSyllIndex).FirstSyllStartTime - PreMotorLag);
    FirstSyll(FirstSyllIndex).FirstSyllDur = FirstSyll(FirstSyllIndex).FirstSyllEndTime - FirstSyll(FirstSyllIndex).FirstSyllStartTime;
    FirstSyll(FirstSyllIndex).FirstSyllAmp = Neural_INR.BoutDetails(BoutNo).Feats(Neural_INR.WithinBoutINs{i}(end)+1, 2);
    FirstSyll(FirstSyllIndex).FirstSyllEnt = Neural_INR.BoutDetails(BoutNo).Feats(Neural_INR.WithinBoutINs{i}(end)+1, 3);
    FirstSyll(FirstSyllIndex).FirstSyllFreq = Neural_INR.BoutDetails(BoutNo).Feats(Neural_INR.WithinBoutINs{i}(end)+1, 4);

    TempIFRIndices =  find((BoutIFR(1,:) >= (FirstSyll(FirstSyllIndex).FirstSyllStartTime - PreMotorLag)) & (BoutIFR(1,:) <= (FirstSyll(FirstSyllIndex).FirstSyllEndTime - PreMotorLag)));
    FirstSyll(FirstSyllIndex).FirstSyllIFR(1,:) = (FirstSyll(FirstSyllIndex).FirstSyllStartTime:1/IFRFs:FirstSyll(FirstSyllIndex).FirstSyllEndTime) - PreMotorLag;
    FirstSyll(FirstSyllIndex).FirstSyllIFR(2,:) = interp1(BoutIFR(1, TempIFRIndices), BoutIFR(2, TempIFRIndices), FirstSyll(FirstSyllIndex).FirstSyllIFR(1,:));
    FirstSyll(FirstSyllIndex).FirstSyllIFR(1,:) = FirstSyll(FirstSyllIndex).FirstSyllIFR(1,:) - (FirstSyll(FirstSyllIndex).FirstSyllStartTime - PreMotorLag);
end

INDetails.INs = IN;
INDetails.FirstSylls = FirstSyll;

disp('Finished Analysis');