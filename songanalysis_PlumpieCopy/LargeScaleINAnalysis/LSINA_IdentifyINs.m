function [INResults] = LSINA_IdentifyINs(BirdParameters)

% This program will identify INs and tabulate their position counted from
% the first or counted from the last before the motif syllable. It will
% also keep track of all INs in between motifs too and give them position
% labels too.

% There are two kinds of ways I've analyzed this:
% 1) The first motif syllable in a bout is considered and what precedes it
% is considered for INs/Syllables that precede motifs
% 2) A user specified common motif is taken and all INs/Syllables before
% these common motifs are considered

% For case 1, I've only looked at the first motif in a bout, whereas for
% case 2, I've systematically looked at all subsequent motifs in a bout too

% =========================================================================
% First a couple of variables that I should initialize

AllBoutStartSyllables = zeros(4, length(BirdParameters.AllOnsets)); % This will have the same length as AllOnsets and will keep track of the syllables at the start of the bout (Case 1) with position too
AllBoutStartINs = zeros(4, length(BirdParameters.AllOnsets)); % Same thing as above for INs

AllCommonMotifPreSyllables = zeros(5, length(BirdParameters.AllOnsets)); % Case 2 - this is all syllables before common motifs with 2 rows for position from first and last and a row indicating whether it is within bout or bout beginning an
AllCommonMotifPreINs = zeros(5, length(BirdParameters.AllOnsets)); % Same as above for INs

AllCommonMotifPositions = zeros(2, length(BirdParameters.AllOnsets)); % This is to keep track of position of a motif within a bout

BoutStarts = find(BirdParameters.BoutLabels == 'Q');
BoutEnds = find(BirdParameters.BoutLabels == 'q');

TotalSyllCount = 0;

for i = 1:length(BoutStarts),
    IndividualBoutLabels{i} = BirdParameters.BoutLabels(BoutStarts(i)+1:BoutEnds(i)-1);
    
    % Case 1: All INs/syllables just before the first motif syllable of a
    % bout
    
    Motifs{i} = zeros(size(IndividualBoutLabels{i}));
    INs{i} = zeros(size(IndividualBoutLabels{i}));
    BoutSyllables{i} = zeros(size(IndividualBoutLabels{i}));
    Others{i} = zeros(size(IndividualBoutLabels{i}));
    
    for j = 1:length(BirdParameters.MotifLabels),
        Motifs{i}(find(IndividualBoutLabels{i} == BirdParameters.MotifLabels(j))) = 1;
    end
    
    for j = 1:length(BirdParameters.INLabels),
        INs{i}(find(IndividualBoutLabels{i} == BirdParameters.INLabels(j))) = 1;
    end
    
    BoutSyllables{i} = ~Motifs{i};
    Others{i} = ~(Motifs{i} + INs{i});
    
    MotifSyllOnsetsOffsets = conv(Motifs{i}, [1 -1], 'same'); 
    MotifOnsets{i} = find(MotifSyllOnsetsOffsets > 0) + 1;
    MotifOffsets{i} = find(MotifSyllOnsetsOffsets < 0);
    
    % Now find INs and syllables before the onset of the first motif with
    % the condition that they are contiguous (for INs) and the largest
    % interval between these contiguous syllables/INs is < the
    % user-specified interval
    
    % First syllables
    Syllables = 1:1:(MotifOnsets{i}(1) - 1);
    TotalNumBoutStartSyllables(i) = length(Syllables);
    if (length(Syllables) > 1)
        Gaps = BirdParameters.AllOnsets(Syllables(2:end) + TotalSyllCount) - BirdParameters.AllOffsets(Syllables(1:end-1) + TotalSyllCount);
        LongGaps = find(Gaps >= BirdParameters.InterINinterval);
        if (isempty(LongGaps))
            BoutStartSyllableIndices{i} = Syllables;
            NumBoutStartSyllables(i) = length(BoutStartSyllableIndices{i});
        else
            BoutStartSyllableIndices{i} = Syllables(LongGaps(end)+1:end);
            NumBoutStartSyllables(i) = length(BoutStartSyllableIndices{i});
        end
    else
        BoutStartSyllableIndices{i} = Syllables;
        NumBoutStartSyllables(i) = length(BoutStartSyllableIndices{i});
    end
       
    % Next INs
    if (~isempty(BoutStartSyllableIndices{i}))
        BoutStartINIndices{i} = find(INs{i}(BoutStartSyllableIndices{i}) > 0) + BoutStartSyllableIndices{i}(1) - 1;
    else
        BoutStartINIndices{i} = [];
    end
    NumBoutStartINs(i) = length(BoutStartINIndices{i});
    if (NumBoutStartINs > 0)
        if (BoutStartINIndices{i}(1) == 1)
            OtherSyllGaps(i) = 1500;
        else
            OtherSyllGaps(i) = BirdParameters.AllOnsets(BoutStartINIndices{i}(1) + TotalSyllCount) - BirdParameters.AllOffsets(BoutStartINIndices{i}(1) - 1 + TotalSyllCount);
        end
    else
        OtherSyllGaps(i) = BirdParameters.AllOnsets(Syllables(end) + 1 + TotalSyllCount) - BirdParameters.AllOffsets(Syllables(end) + TotalSyllCount);
    end
    if (NumBoutStartSyllables(i) > 0)
        AllBoutStartSyllables(1, BoutStartSyllableIndices{i} + TotalSyllCount) = 1;
        AllBoutStartSyllables(2, BoutStartSyllableIndices{i} + TotalSyllCount) = i;
        AllBoutStartSyllables(3, BoutStartSyllableIndices{i} + TotalSyllCount) = 1:1:length(BoutStartSyllableIndices{i});
        AllBoutStartSyllables(4, BoutStartSyllableIndices{i} + TotalSyllCount) = -length(BoutStartSyllableIndices{i}):1:-1;
    end
    
    if (NumBoutStartINs(i) > 0)
        AllBoutStartINs(1, BoutStartINIndices{i} + TotalSyllCount) = 1;
        AllBoutStartINs(2, BoutStartINIndices{i} + TotalSyllCount) = i;
        AllBoutStartINs(3, BoutStartINIndices{i} + TotalSyllCount) = 1:1:length(BoutStartINIndices{i});
        AllBoutStartINs(4, BoutStartINIndices{i} + TotalSyllCount) = -length(BoutStartINIndices{i}):1:-1;
    end
    
    % Case 2: All common motifs in a bout and the INs/syllables before that
    % In case there is more than one common motif, then I should check to
    % see if one of the common motifs is a subset of the other
    
    if (length(BirdParameters.CommonMotifs) == 1)
        CommonMotifOnsets = strfind(IndividualBoutLabels{i}, BirdParameters.CommonMotifs{1});
        CommonMotifOffsets = CommonMotifOnsets + length(BirdParameters.CommonMotifs{1}) - 1;
    else
        CommonMotifOnsets = [];
        CommonMotifOffsets = [];
        for j = 1:length(BirdParameters.CommonMotifs),
            TempMotifOnsets = strfind(IndividualBoutLabels{i}, BirdParameters.CommonMotifs{j});
            CommonMotifOnsets = [CommonMotifOnsets TempMotifOnsets(:)'];
            CommonMotifOffsets = [CommonMotifOffsets (TempMotifOnsets(:)' + length(BirdParameters.CommonMotifs{j}) - 1)];
        end
        [CommonMotifOnsets, SortedIndices] = sort(CommonMotifOnsets);
        CommonMotifOffsets = CommonMotifOffsets(SortedIndices);
        
        [CommonMotifOffsets, UniqueIndices] = unique(CommonMotifOffsets);
        CommonMotifOnsets = CommonMotifOnsets(UniqueIndices);
    end     
    for j = 1:length(CommonMotifOnsets),
        AllCommonMotifPositions(1,(CommonMotifOnsets(j):CommonMotifOffsets(j)) + TotalSyllCount) = j;
        AllCommonMotifPositions(2,(CommonMotifOnsets(j):CommonMotifOffsets(j)) + TotalSyllCount) = i;
        if (j == 1)
            Syllables = 1:1:CommonMotifOnsets(1)-1;
        else
            Syllables = (CommonMotifOffsets(j-1)+1):1:(CommonMotifOnsets(j)-1);
        end
        if (length(Syllables) > 1)
            Gaps = BirdParameters.AllOnsets(Syllables(2:end) + TotalSyllCount) - BirdParameters.AllOffsets(Syllables(1:end-1) + TotalSyllCount);
            LongGaps = find(Gaps >= BirdParameters.InterINinterval);
            if (isempty(LongGaps))
                CommonMotifStartSyllableIndices{i}{j} = Syllables;
                NumCommonMotifStartSyllables{i}(j) = length(CommonMotifStartSyllableIndices{i}{j});
            else
                CommonMotifStartSyllableIndices{i}{j} = Syllables(LongGaps(end)+1:end);
                NumCommonMotifStartSyllables{i}(j) = length(CommonMotifStartSyllableIndices{i}{j});
            end
        else
            CommonMotifStartSyllableIndices{i}{j} = Syllables;
            NumCommonMotifStartSyllables{i}(j) = length(CommonMotifStartSyllableIndices{i}{j});
        end
        
        if (~isempty(CommonMotifStartSyllableIndices{i}{j}))
            CommonMotifStartINIndices{i}{j} = find(INs{i}(CommonMotifStartSyllableIndices{i}{j}) > 0) + CommonMotifStartSyllableIndices{i}{j}(1) - 1;
        else
            CommonMotifStartINIndices{i}{j} = [];
        end
        NumCommonMotifStartINs{i}(j) = length(CommonMotifStartINIndices{i}{j});
    end
    
    for j = 1:length(NumCommonMotifStartSyllables{i}),
        if (NumCommonMotifStartSyllables{i}(j) > 0)    
            AllCommonMotifPreSyllables(1, CommonMotifStartSyllableIndices{i}{j} + TotalSyllCount) = 1;
            AllCommonMotifPreSyllables(2, CommonMotifStartSyllableIndices{i}{j} + TotalSyllCount) = i;
            AllCommonMotifPreSyllables(3, CommonMotifStartSyllableIndices{i}{j} + TotalSyllCount) = j;
            AllCommonMotifPreSyllables(4, CommonMotifStartSyllableIndices{i}{j} + TotalSyllCount) = 1:1:length(CommonMotifStartSyllableIndices{i}{j});
            AllCommonMotifPreSyllables(5, CommonMotifStartSyllableIndices{i}{j} + TotalSyllCount) = -length(CommonMotifStartSyllableIndices{i}{j}):1:-1;
        end
    end

    for j = 1:length(NumCommonMotifStartINs{i}),
        if (NumCommonMotifStartINs{i}(j) > 0)    
            AllCommonMotifPreINs(1, CommonMotifStartINIndices{i}{j} + TotalSyllCount) = 1;
            AllCommonMotifPreINs(2, CommonMotifStartINIndices{i}{j} + TotalSyllCount) = i;
            AllCommonMotifPreINs(3, CommonMotifStartINIndices{i}{j} + TotalSyllCount) = j;
            AllCommonMotifPreINs(4, CommonMotifStartINIndices{i}{j} + TotalSyllCount) = 1:1:length(CommonMotifStartINIndices{i}{j});
            AllCommonMotifPreINs(5, CommonMotifStartINIndices{i}{j} + TotalSyllCount) = -length(CommonMotifStartINIndices{i}{j}):1:-1;
        end
    end
    TotalSyllCount = TotalSyllCount + length(IndividualBoutLabels{i});
end

INResults.IndividualBoutLabels = IndividualBoutLabels;
INResults.Motifs = Motifs;
INResults.INs = INs;
INResults.Syllables = BoutSyllables;
INResults.Others = Others;
    
INResults.MotifOnsets = MotifOnsets;
INResults.MotifOffsets = MotifOffsets;
INResults.BoutStartSyllableIndices = BoutStartSyllableIndices;
INResults.NumBoutStartSyllables = NumBoutStartSyllables;
INResults.BoutStartINIndices = BoutStartINIndices;
INResults.NumBoutStartINs = NumBoutStartINs;
INResults.AllBoutStartSyllables = AllBoutStartSyllables;
INResults.AllBoutStartINs = AllBoutStartINs;
INResults.TotalNumBoutStartSyllables = TotalNumBoutStartSyllables;

INResults.Gaps = OtherSyllGaps;

INResults.CommonMotifStartSyllableIndices = CommonMotifStartSyllableIndices;
INResults.NumCommonMotifStartSyllables = NumCommonMotifStartSyllables;
INResults.CommonMotifStartINIndices = CommonMotifStartINIndices;
INResults.NumCommonMotifStartINs = NumCommonMotifStartINs;

INResults.AllCommonMotifPreINs = AllCommonMotifPreINs;
INResults.AllCommonMotifPreSyllables = AllCommonMotifPreSyllables;
INResults.AllCommonMotifPositions = AllCommonMotifPositions;

disp('Finished');