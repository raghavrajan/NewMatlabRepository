function [IN, Motif, Bout] = CSIdentifyINs(Labels, MotifSylls, INLabels)

LabelSymbols = Labels;
LabelSymbols(regexp(Labels, eval(['''[', MotifSylls, ']''']))) = 'M';
LabelSymbols(regexp(Labels, eval(['''[', INLabels, ']''']))) = 'I';

MotifLabels = (LabelSymbols == 'M');
% MotifLabels is now a logical variable with 1's and 0's corresponding to
% locations of Ms and non-Ms. This can be converted to a double by doing
% +MotifLabels - any arithmetic operation


MotifTrans = conv(+MotifLabels, [1 -1]);
Motif.Starts = find(MotifTrans == 1);
Motif.Ends = find(MotifTrans == -1) - 1;
Motif.Indices = find(+MotifLabels == 1);

Bout.Starts = find(Labels == 'Q');
Bout.Starts = Bout.Starts(:);

Bout.Ends = find(Labels == 'q');
Bout.Ends = Bout.Ends(:);

% Now each motif has to be assigned to a bout and within that bout, the
% motif has to be given a position as the 1st or the 2nd or ..... motif in
% the bout.

Motif.IndicesLocations = repmat(Motif.Indices, length(Bout.Starts), 1);
Motif.IndicesLocations = (Motif.IndicesLocations - repmat(Bout.Starts, 1, size(Motif.IndicesLocations, 2))).*(repmat(Bout.Ends, 1, size(Motif.IndicesLocations, 2)) - Motif.IndicesLocations);
Motif.IndicesLocations = +(Motif.IndicesLocations > 0);
Motif.IndicesLocations = (Motif.IndicesLocations).*(repmat((1:1:length(Bout.Starts))', 1, size(Motif.IndicesLocations,2)));
Motif.IndicesLocations = sum(Motif.IndicesLocations);

Motif.StartsLocations = repmat(Motif.Starts, length(Bout.Starts), 1);
Motif.StartsLocations = (Motif.StartsLocations - repmat(Bout.Starts, 1, size(Motif.StartsLocations, 2))).*(repmat(Bout.Ends, 1, size(Motif.StartsLocations, 2)) - Motif.StartsLocations);
Motif.StartsLocations = +(Motif.StartsLocations > 0);
Motif.StartsLocations = (Motif.StartsLocations).*(repmat((1:1:length(Bout.Starts))', 1, size(Motif.StartsLocations,2)));
Motif.StartsLocations = sum(Motif.StartsLocations);
Motif.EndsLocations = Motif.StartsLocations;

Motif.BoutBeginningMotifs = find(diff(Motif.StartsLocations));
Motif.BoutBeginningMotifs = [1 (Motif.BoutBeginningMotifs(:)' + 1)];
Motif.WithinBoutMotifs = setdiff((1:1:length(Motif.Starts)), Motif.BoutBeginningMotifs);

Motif.MotifNumber = 1:1:length(Motif.Starts);

INSequences = (LabelSymbols == 'I');
INTrans = conv(+INSequences, [1 -1]);
IN.Starts = find(INTrans == 1);
IN.Ends = find(INTrans == -1) - 1;
IN.Indices = find(+INSequences == 1);

% Now each IN has to be assigned to a bout and within that bout, the
% IN has to be given a position as the 1st or the 2nd or ..... IN in
% the bout. The ordering has to be done both assuming the first IN is
% common and assuming the last IN is common.

IN.IndicesLocations = repmat(IN.Indices, length(Bout.Starts), 1);
IN.IndicesLocations = (IN.IndicesLocations - repmat(Bout.Starts, 1, size(IN.IndicesLocations, 2))).*(repmat(Bout.Ends, 1, size(IN.IndicesLocations, 2)) - IN.IndicesLocations);
IN.IndicesLocations = +(IN.IndicesLocations > 0);
IN.IndicesLocations = (IN.IndicesLocations).*(repmat((1:1:length(Bout.Starts))', 1, size(IN.IndicesLocations,2)));
IN.IndicesLocations = sum(IN.IndicesLocations);

IN.StartsLocations = repmat(IN.Starts, length(Bout.Starts), 1);
IN.StartsLocations = (IN.StartsLocations - repmat(Bout.Starts, 1, size(IN.StartsLocations, 2))).*(repmat(Bout.Ends, 1, size(IN.StartsLocations, 2)) - IN.StartsLocations);
IN.StartsLocations = +(IN.StartsLocations > 0);
IN.StartsLocations = (IN.StartsLocations).*(repmat((1:1:length(Bout.Starts))', 1, size(IN.StartsLocations,2)));
IN.StartsLocations = sum(IN.StartsLocations);

IN.MotifLocations = repmat(IN.Indices, length(Motif.Starts), 1);
IN.MotifLocations = (repmat(Motif.Starts', 1, size(IN.MotifLocations, 2)) - IN.MotifLocations).*(IN.MotifLocations - repmat([1 Motif.Ends(1:end-1)]', 1, size(IN.MotifLocations, 2)));
IN.MotifLocations = +(IN.MotifLocations > 0);
IN.MotifLocations = (IN.MotifLocations).*(repmat((1:1:length(Motif.Starts))', 1, size(IN.MotifLocations,2)));

for i = 1:length(Motif.Starts),
    TempLocations = find(IN.MotifLocations(i,:));
    Temp = diff([IN.Indices(TempLocations) Motif.Starts(i)]);
    Temp = find(Temp > 1);
    if (~isempty(Temp))
        IN.MotifLocations(i, TempLocations(1:Temp(end))) = 0;
    end
    IN.NumINs(i) = length(find(IN.MotifLocations(i,:)));
 
    IN.PosFromFirst(i,:) = IN.MotifLocations(i,:);
    IN.PosFromFirst(i,find(IN.PosFromFirst(i,:))) = 1:1:IN.NumINs(i);
    
    IN.PosFromLast(i,:) = IN.MotifLocations(i,:);
    IN.PosFromLast(i,find(IN.PosFromLast(i,:))) = -IN.NumINs(i):1:-1;
end

disp('Finished identifying INs');