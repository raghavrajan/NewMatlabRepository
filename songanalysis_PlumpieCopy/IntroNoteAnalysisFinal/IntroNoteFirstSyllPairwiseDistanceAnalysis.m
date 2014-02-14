function [R, P] = IntroNoteFirstSyllPairwiseDistanceAnalysis(IntroNoteResults, Motif, varargin)

Index = 0;
LastINIndex = 0;
for i = 1:length(IntroNoteResults.NoofINs),
    Index = Index + 1;
    if (IntroNoteResults.NoofINs(i) > 0)
        INs = IntroNoteResults.INs{i};
        LastINIndex = LastINIndex + 1;
        LastINFeatures(LastINIndex,:) = IntroNoteResults.BoutDetails(i).Feats(INs(end),1:4);
        FirstSyllFeatures(LastINIndex,:) = IntroNoteResults.BoutDetails(i).Feats(IntroNoteResults.MotifStartIndex(i),1:4);
        FirstSyllBoutPos(LastINIndex) = 0;
    end
    AllFirstSyllFeatures(Index,:) = IntroNoteResults.BoutDetails(i).Feats(IntroNoteResults.MotifStartIndex(i),1:4);
    AllFirstSyllBoutPos(Index) = 0;
    NoofINs(Index, :) = IntroNoteResults.NoofINs(i);
    
%     Matches = strfind(IntroNoteResults.BoutDetails(i).labels, Motif);
%     if (~isempty(find(Matches == IntroNoteResults.MotifStartIndex(i))))
%        MotifOnset = IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.MotifStartIndex(i));
%         MotifOffset = IntroNoteResults.BoutDetails(i).offsets(IntroNoteResults.MotifStartIndex(i) + length(Motif) - 1);
%         MotifWaveform = TempSong(find((Time >= (MotifOnset - 0.005)) & (Time <= (MotifOffset + 0.005))));
%         cd(OutputDir);
%         wavwrite(MotifWaveform, Fs, 16, [OutputFileName, '.Motif.', num2str(Index), '.wav']);
%     end
end

for i = 1:size(IntroNoteResults.WithinBoutNoofINs,1),
    Index = Index + 1;
    if (IntroNoteResults.WithinBoutNoofINs(i,1) > 0)
        INs = IntroNoteResults.WithinBoutINs{i};
        LastINIndex = LastINIndex + 1;
        LastINFeatures(LastINIndex,:) = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).Feats(INs(end),1:4);
        FirstSyllFeatures(LastINIndex,:) = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).Feats(IntroNoteResults.WithinBoutNoofINs(i, 4),1:4);
        FirstSyllBoutPos(LastINIndex) = 1;
    end
    AllFirstSyllFeatures(Index,:) = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).Feats(IntroNoteResults.WithinBoutNoofINs(i, 4),1:4);
    AllFirstSyllBoutPos(Index) = 1;
    NoofINs(Index, :) = IntroNoteResults.WithinBoutNoofINs(i,1);
%     Matches = strfind(IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).labels, Motif);
%     if (~isempty(find(Matches == IntroNoteResults.WithinBoutNoofINs(i, 4))))
%         MotifOnset = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).onsets(IntroNoteResults.WithinBoutNoofINs(i, 4));
%         MotifOffset = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).offsets(IntroNoteResults.WithinBoutNoofINs(i, 4) + length(Motif) - 1);
%         MotifWaveform = TempSong(find((Time >= (MotifOnset - 0.005)) & (Time <= (MotifOffset + 0.005))));
%         cd(OutputDir);
%         wavwrite(MotifWaveform, Fs, 16, [OutputFileName, '.Motif.', num2str(Index), '.wav']);
%     end
end

[NaNLastINs,j] = find(isnan(LastINFeatures));
if (~isempty(NaNLastINs))
    LastINFeatures(NaNLastINs, :) = [];
    FirstSyllFeatures(NaNLastINs, :) = [];
    FirstSyllBoutPos(NaNLastINs) = [];
end

[NaNFirstSylls, j] = find(isnan(FirstSyllFeatures));
if (~isempty(NaNFirstSylls))
    LastINFeatures(NaNFirstSylls, :) = [];
    FirstSyllFeatures(NaNFirstSylls, :) = [];
    FirstSyllBoutPos(NaNFirstSylls) = [];
end

[NaNFirstSylls, j] = find(isnan(AllFirstSyllFeatures));
if (~isempty(NaNFirstSylls))
    NoofINs(NaNFirstSylls, :) = [];
    AllFirstSyllFeatures(NaNFirstSylls, :) = [];
    AllFirstSyllBoutPos(NaNFirstSylls) = [];
end

PairWiseLastINDist = pdist2(LastINFeatures, LastINFeatures, 'mahalanobis', nancov(LastINFeatures));
PairWiseLastINDist = PairWiseLastINDist(find(tril(PairWiseLastINDist, -1)));
 
PairWiseFirstSyllDist = pdist2(FirstSyllFeatures, FirstSyllFeatures, 'mahalanobis', nancov(AllFirstSyllFeatures));
PairWiseFirstSyllDist = PairWiseFirstSyllDist(find(tril(PairWiseFirstSyllDist, -1)));

PairWiseWithinBoutFirstSyllDist = pdist2(FirstSyllFeatures(find(FirstSyllBoutPos == 1),:), FirstSyllFeatures, 'mahalanobis', nancov(AllFirstSyllFeatures));
PairWiseWithinBoutFirstSyllDist = PairWiseWithinBoutFirstSyllDist(find(tril(PairWiseWithinBoutFirstSyllDist, -1)));

NoINPairWiseFirstSyllDist = pdist2(AllFirstSyllFeatures(find(NoofINs == 0), :), FirstSyllFeatures, 'mahalanobis', nancov(AllFirstSyllFeatures));
NoINPairWiseFirstSyllDist = NoINPairWiseFirstSyllDist(find(tril(NoINPairWiseFirstSyllDist, -1)));

disp(['Mean distance between first syllables with INs: ', num2str(mean(PairWiseFirstSyllDist)), ' and between first syllables with INs and without INs:', num2str(mean(NoINPairWiseFirstSyllDist))]);

%figure;
%[p, anovatab, stats] = kruskalwallis([PairWiseFirstSyllDist; NoINPairWiseFirstSyllDist], [ones(size(PairWiseFirstSyllDist))*1; ones(size(NoINPairWiseFirstSyllDist))*2]);

figure
plot(FirstSyllFeatures(find(FirstSyllBoutPos == 1),2), FirstSyllFeatures(find(FirstSyllBoutPos == 1),3), 'ro');
hold on;
plot(AllFirstSyllFeatures(find(NoofINs == 0),2), AllFirstSyllFeatures(find(NoofINs == 0),3), 'bo');
xlabel('Log amplitude');
ylabel('Wiener entropy');

[r, p] = corrcoef(PairWiseLastINDist, PairWiseFirstSyllDist);
R = r(1,2);
P = p(1,2);

figure;
plot(PairWiseLastINDist, PairWiseFirstSyllDist, 'k+', 'MarkerSize', 3);
axis tight;
temp = axis;
hold on;
temp = [min([temp(1) temp(3)]) max([temp(2) temp(4)]) min([temp(1) temp(3)]) max([temp(2) temp(4)])];
axis(temp);
plot([temp(1) temp(2)], [temp(1) temp(2)], 'r--', 'LineWidth', 2);

disp('Finished feature analysis');
