function [BoutStartSyll, FirstMotifSyll, PreFirstMotifSyll, FirstMotifSyllIndex] = LSINA_CalculateBoutStartingSyllStats(BirdParameters)

BoutStartIndices = find(BirdParameters.BoutLabels == 'Q');

% First get the bout starting syll
BoutStartSyll = BirdParameters.BoutLabels(BoutStartIndices + 1);

% Now find the first motif syll
for i = 1:length(BoutStartIndices),
    Flag = 1;
    Index = BoutStartIndices(i);
    while (Flag == 1)
        Index = Index + 1;
        if (~isempty(find(BirdParameters.MotifLabels == BirdParameters.BoutLabels(Index))))
            Flag = 0;
            FirstMotifSyllIndex(i) = Index;
            FirstMotifSyll(i) = BirdParameters.BoutLabels(Index);
            PreFirstMotifSyll(i) = BirdParameters.BoutLabels(Index - 1);
        end
    end
end