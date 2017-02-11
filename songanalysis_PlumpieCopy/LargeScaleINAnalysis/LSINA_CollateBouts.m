function [BoutLabels] = LSINA_CollateBouts(BirdParameters)

% Now that bouts have been identified, collate all bouts together for
% calculating transition probabilities
% Basically, add a 'Q' at the beginning of each bout and a 'q' at the end
% of each bout

BoutLabels = [];
if (BirdParameters.Continuousdata == 0)
    for j = 1:length(BirdParameters.Bouts),
        if (BirdParameters.Bouts(j) == 1)
            for k = 1:size(BirdParameters.BoutIndices{j},1),
                BoutLabels = [BoutLabels 'Q' BirdParameters.NoteInfo{j}.labels(BirdParameters.BoutIndices{j}(k,1):BirdParameters.BoutIndices{j}(k,2)), 'q'];
            end
        end
    end
else
    for j = 1:size(BirdParameters.BoutIndices,1),
        BoutLabels = [BoutLabels 'Q' BirdParameters.AllLabels(BirdParameters.BoutIndices(j,1):BirdParameters.BoutIndices(j,2)), 'q'];
    end
end
