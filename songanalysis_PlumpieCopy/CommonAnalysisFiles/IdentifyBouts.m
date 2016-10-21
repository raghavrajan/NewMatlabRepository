function [Bouts] = IdentifyBouts(BirdParameters)

% Now based on whether the data is continuous or not, identify bouts of
% song using SyllableData, the field of BirdParameters, that has all the
% data about all the syllables
Bouts = [];
if (BirdParameters.Continuousdata == 0)
    UniqueFileIndices = unique(BirdParameters.SyllableData(:,2));
    for i = 1:length(UniqueFileIndices),
        % Within each file, I should check if there are gaps longer than
        % the specified inter-bout interval. Then do as follows:
        % a) if there are no long gaps, then check if there is enough space
        % at the beginning and end of the bout
        % b) if there are long gaps, then accordingly sort out bouts and
        % check if there is enough space at the beginning and end of each
        % bout
        % Also check that each bout has motif syllables in it too for it to
        % be classified as a bout
        
        FileIndices = find(BirdParameters.SyllableData(:,2) == UniqueFileIndices(i));
        Gaps = BirdParameters.SyllableData(FileIndices(2:end),4) - BirdParameters.SyllableData(FileIndices(1:end-1),5);
        LongGaps = find(Gaps >= BirdParameters.Interboutinterval);
        if (isempty(LongGaps))
            if (~isempty(strfind(char(BirdParameters.SyllableData(FileIndices,1))', BirdParameters.CommonMotifs{1})))
                if ((BirdParameters.SyllableData(FileIndices(1),4) >= BirdParameters.Interboutinterval) && (BirdParameters.SyllableData(FileIndices(end),5) <= (BirdParameters.FileLen(i) - BirdParameters.Interboutinterval)))
                    Bouts(end+1,:) = [FileIndices(1) FileIndices(end)];
                end
            end
        else
            for k = 1:length(LongGaps),
                if (k == 1)
                    if (~isempty(strfind(char(BirdParameters.SyllableData(FileIndices(1:LongGaps(k)),1))', BirdParameters.CommonMotifs{1})))
                        if (BirdParameters.SyllableData(FileIndices(1),4) >= BirdParameters.Interboutinterval)
                            Bouts(end+1,:) = [FileIndices(1) FileIndices(LongGaps(k))];
                        end
                    end
                else
                    if (~isempty(strfind(char(BirdParameters.SyllableData(FileIndices(LongGaps(k-1)+1:LongGaps(k)),1))', BirdParameters.CommonMotifs{1})))
                        Bouts(end+1,:) = [FileIndices(LongGaps(k-1)+1) FileIndices(LongGaps(k))];
                    end
                end
            end
            if (~isempty(strfind(char(BirdParameters.SyllableData(FileIndices(LongGaps(end)+1:end),1))', BirdParameters.CommonMotifs{1})))
                if (BirdParameters.SyllableData(FileIndices(end),5) <= (BirdParameters.FileLen(i) - BirdParameters.Interboutinterval))
                    Bouts(end+1,:) = [FileIndices(LongGaps(end) + 1) FileIndices(end)];
                end
            end
        end
    end
else
    Gaps = BirdParameters.SyllableData(2:end, 6) - BirdParameters.SyllableData(1:end-1,7);
    LongGaps = find(Gaps >= BirdParameters.Interboutinterval);

    if (isempty(LongGaps))
        if (~isempty(strfind(char(BirdParameters.SyllableData(:,1))', BirdParameters.CommonMotifs{1})))
            if ((BirdParameters.SyllableData(1,6) >= BirdParameters.Interboutinterval) && (BirdParameters.SyllableData(end,7) <= (sum(BirdParameters.FileLen) - BirdParameters.Interboutinterval)))
                Bouts(end+1,:) = [1 size(BirdParameters.SyllableData,1)];
            end
        end
    else
        for k = 1:length(LongGaps),
            if (k == 1)
                if (~isempty(strfind(char(BirdParameters.SyllableData(1:LongGaps(k),1))', BirdParameters.CommonMotifs{1})))
                    if (BirdParameters.SyllableData(1,6) >= BirdParameters.Interboutinterval)
                        Bouts(end+1,:) = [1 LongGaps(k)];
                    end
                end
            else
                if (~isempty(strfind(char(BirdParameters.SyllableData((LongGaps(k-1)+1):LongGaps(k),1))', BirdParameters.CommonMotifs{1})))
                    Bouts(end+1,:) = [(LongGaps(k-1)+1) LongGaps(k)];
                end
            end
        end
        if (~isempty(strfind(char(BirdParameters.SyllableData(LongGaps(end)+1:end,1))', BirdParameters.CommonMotifs{1})))
            if (BirdParameters.SyllableData(end,7) <= (sum(BirdParameters.FileLen) - BirdParameters.Interboutinterval))
                Bouts(end+1,:) = [(LongGaps(end)+1) size(BirdParameters.SyllableData,1)];
            end
        end
    end
end
