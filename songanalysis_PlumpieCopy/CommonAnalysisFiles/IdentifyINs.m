function [INs, Motifs] = IdentifyINs(BirdParameters)

% This program will identify INs and tabulate their position counted from
% the first or counted from the last before the motif syllable. It will
% also keep track of all INs in between motifs too and give them position
% labels too.

% Things I want to keep track of for the matrix Motifs
% a) column 1 - whether it is a motif syllable or not
% b) column 2 - which motif in the bout it is
% 
% Things I want to keep track of for the matrix INs
% a) column 1 - whether it is an IN or not
% b) column 2 - whether it is at the beginning (1) or within bout (2)
% c) column 3 - position of IN relative to first IN
% d) column 4 - position of IN relative to last IN

% Finally one last variable called IN_NonIN_Bout to keep track of whether
% the out has only INs at the beginning (1) or has other types of syllables (0) too
% before the first motif syllable 

% =========================================================================

INs = zeros(size(BirdParameters.SyllableData,1),4);
Motifs = zeros(size(BirdParameters.SyllableData,1),2);
IN_NonIN_Bout = zeros(size(BirdParameters.Bouts,1),1);

for i = 1:size(BirdParameters.Bouts,1),
   BoutLabels = char(BirdParameters.SyllableData(Bouts(i,1):Bouts(i,2)))';
   
   for j = 1:length(BoutLabels),
       if (BirdParameters.MotifLabels == BoutLabels(j))
           Motifs(Bouts(i,1) + j - 1, 1) = 1;
       end
   end
   
   % Now find all motifs in the bout
    if (length(BirdParameters.CommonMotifs) == 1)
       CommonMotifIndices = strfind(BoutLabels, BirdParameters.CommonMotifs{1});
       for k = 1:length(CommonMotifIndicess),
           Motifs(Bouts(i,1) + CommonMotifIndicess(k) - 1:Bouts(i,1) - 1 + length(BirdParameters.CommonMotifs{1}) - 1) = k;
       end
    else
        CommonMotifOnsets = [];
        CommonMotifOffsets = [];
        for k = 1:length(BirdParameters.CommonMotifs),
            TempMotifOnsets = strfind(BoutLabels, BirdParameters.CommonMotifs{k});
            CommonMotifOnsets = [CommonMotifOnsets TempMotifOnsets(:)'];
            CommonMotifOffsets = [CommonMotifOffsets (TempMotifOnsets(:)' + length(BirdParameters.CommonMotifs{j}) - 1)];
        end
        [CommonMotifOnsets, SortedIndices] = sort(CommonMotifOnsets);
        CommonMotifOffsets = CommonMotifOffsets(SortedIndices);

        [CommonMotifOffsets, UniqueIndices] = unique(CommonMotifOffsets);
        CommonMotifOnsets = CommonMotifOnsets(UniqueIndices);
        for k = 1:length(CommonMotifOnsets),
            Motifs(Bouts(i,1) - 1 + CommonMotifOnsets(k):Bouts(i,1) - 1 + CommonMotifOffsets(k)) = k;
        end
    end     
    
    for j = 1:length(BoutLabels),
       if (BirdParameters.INLabels == BoutLabels(j))
           INs(Bouts(i,1) + j - 1, 1) = 1;
           if (j < CommonMotifOnsets(1))
               INs(Bouts(i,1) + j - 1, 2) = 1;
           else
               INs(Bouts(i,1) + j - 1, 2) = 2;
           end
       end
    end

    % Now I have to figure out sequences of INs such that they are
    % continuous and don't have more than 500 ms between them.
    for j = 1:length(CommonMotifOnsets),
        if (j == 1)
            Syllables = 1:(CommonMotifOnsets(j) - 1);
            INSyllables = zeros(size(Syllables));
            for k = 1:length(Syllables),
                if (BirdParameters.INLabels == BoutLabels(Syllables))
                    INSyllables(k) = 1;
                end
            end
            INSyllables = find(INSyllables);
            
            if (BirdParameters.Continuousdata == 0)
                Gaps = BirdParameters.SyllableData(Bouts(i,1) - 1 + Syllables(2):Bouts(i,1) - 1 + Syllables(end), 4) - BirdParameters.SyllableData(Bouts(i,1) - 1 + Syllables(2):Bouts(i,1) - 1 + Syllables(end), 5);
                LongGaps = find(Gaps >= BirdParameters.InterINinterval);
                
                
end
   