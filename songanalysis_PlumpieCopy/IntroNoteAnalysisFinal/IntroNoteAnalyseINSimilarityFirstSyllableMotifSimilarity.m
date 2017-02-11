function [IN_A_Sim, IN_A_OnsetTime, IN_Motif_Sim] = IntroNoteAnalyseINSimilarityFirstSyllableMotifSimilarity(IntroNoteResults, CSVFile, Motif)

if (isempty(CSVFile))
    IN_A_Sim(1,:) = [NaN NaN];
    IN_Motif_Sim(1,:) = [NaN NaN];
    IN_A_OnsetTime(1,:) = [NaN NaN];
    return;
end
Fid = fopen(CSVFile, 'r');

Data = textscan(Fid, '%d%s%s%f%f%f%f%f%f%f%f%f%f%f%f', 'DeLimiter', ',', 'HeaderLines', 1);
Data([4:7 11:14]) = [];
fclose(Fid);

PresentDir = pwd;

Index = 0;
for i = 1:length(IntroNoteResults.NoofINs),
    Index = Index + 1;
    if (IntroNoteResults.NoofINs(i) > 0)
        INs = IntroNoteResults.INs{i};
        INOnset = IntroNoteResults.BoutDetails(i).onsets(INs(end));
        INOffset = IntroNoteResults.BoutDetails(i).offsets(INs(end));
        INDur(Index) = INOffset - INOnset;
        
        FirstSyllOnset = IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.MotifStartIndex(i));
        FirstSyllOffset = IntroNoteResults.BoutDetails(i).offsets(IntroNoteResults.MotifStartIndex(i));
        
        FirstSyllOnsetTime(Index,:) = [(FirstSyllOnset - INOffset) (FirstSyllOnset - INOnset)];
        FirstSyllDur(Index) = FirstSyllOffset - FirstSyllOnset;
        
        Matches = strfind(IntroNoteResults.BoutDetails(i).labels, Motif);
        if (~isempty(find(Matches == IntroNoteResults.MotifStartIndex(i))))
            MotifOnset = IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.MotifStartIndex(i));
            MotifOffset = IntroNoteResults.BoutDetails(i).offsets(IntroNoteResults.MotifStartIndex(i) + length(Motif) - 1);
            
            MotifOnsetTime(Index,:) = [(MotifOnset - INOffset) (MotifOnset - INOnset)];
            MotifDur(Index) = MotifOffset - MotifOnset;
        end
    end
end

for i = 1:size(IntroNoteResults.WithinBoutNoofINs,1),
    Index = Index + 1;
    if (IntroNoteResults.WithinBoutNoofINs(i,1) > 0)
        INs = IntroNoteResults.WithinBoutINs{i};
        INOnset = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).onsets(INs(end));
        INOffset = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).offsets(INs(end));
        INDur(Index) = INOffset - INOnset;
     
        FirstSyllOnset = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).onsets(IntroNoteResults.WithinBoutNoofINs(i, 4));
        FirstSyllOffset = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).offsets(IntroNoteResults.WithinBoutNoofINs(i, 4));
        
        FirstSyllOnsetTime(Index,:) = [(FirstSyllOnset - INOffset) (FirstSyllOnset - INOnset)];
        FirstSyllDur(Index) = FirstSyllOffset - FirstSyllOnset;
        
        Matches = strfind(IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).labels, Motif);
        if (~isempty(find(Matches == IntroNoteResults.WithinBoutNoofINs(i, 4))))
            MotifOnset = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).onsets(IntroNoteResults.WithinBoutNoofINs(i, 4));
            MotifOffset = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).offsets(IntroNoteResults.WithinBoutNoofINs(i, 4) + length(Motif) - 1);

            MotifOnsetTime(Index,:) = [(MotifOnset - INOffset) (MotifOnset - INOnset)];
            MotifDur(Index) = MotifOffset - MotifOnset;
        end
    end
end

As = find(cellfun(@length, strfind(Data{3}, '.a.')));
for i = 1:length(As),
    StartIndex = strfind(Data{3}{As(i)}, '.a.');
    StartIndex = StartIndex + 3;
    EndIndex = strfind(Data{3}{As(i)}, '.wav');
    EndIndex = EndIndex - 1;
    AIndex(i) = str2double(Data{3}{As(i)}(StartIndex:EndIndex));
end

for i = 1:max(AIndex),
    TempIN = find(cellfun(@length, strfind(Data{3}, ['.IN.', num2str(i), '.wav'])));
    if (~isempty(TempIN))
        Similarity(i,1) = Data{4}(TempIN(1));
    else
        Similarity(i,1) = NaN;
    end
    TempIN = find(cellfun(@length, strfind(Data{3}, ['.a.', num2str(i), '.wav'])));
    if (~isempty(TempIN))
        Similarity(i,2) = Data{4}(TempIN(1));
    else
        Similarity(i,2) = NaN;    
    end
    TempIN = find(cellfun(@length, strfind(Data{3}, ['.Motif.', num2str(i), '.wav'])));
    if (~isempty(TempIN))
        Similarity(i,3) = Data{4}(TempIN(1));
    else
        Similarity(i,3) = NaN;
    end
end

Similarity(1,:) = [];

ZeroINSim = find(Similarity(:,1) == 0);
Similarity(ZeroINSim,:) = [];

ZeroASim = find(Similarity(:,2) == 0);
Similarity(ZeroASim,:) = [];

ZeroMotifSim = find(Similarity(:,3) == 0);
Similarity(ZeroMotifSim,:) = [];

INs_and_As = find((~isnan(Similarity(:,1))) & (~isnan(Similarity(:,2))));
INs_and_Motifs = find((~isnan(Similarity(:,1))) & (~isnan(Similarity(:,3))));

figure;
subplot(2,1,1);
hold on;
plot(Similarity(INs_and_As,1), Similarity(INs_and_As,2), 'ko');
axis tight;
temp = axis;
axis([temp(1:3) temp(4)+3]);

NoINs = find(isnan(Similarity(:,1)));
if (~isempty(NoINs))
    plot(ones(size(NoINs))*(temp(4)+1.5), Similarity(NoINs, 2), 'ro');
end
axis tight;

subplot(2,1,2);
plot(Similarity(INs_and_Motifs, 1), Similarity(INs_and_Motifs, 3), 'ko');
axis tight;
temp = axis;
axis([temp(1:3) temp(4)+3]);

hold on;
if (~isempty(NoINs))
    plot(ones(size(NoINs))*(temp(4)+1.5), Similarity(NoINs, 3), 'ro');
end
axis tight;
[r, p] = corrcoef(Similarity(INs_and_As, 1), Similarity(INs_and_As,2));
disp(['Correlation between IN and a: r = ', num2str(r(1,2)), ' and p = ', num2str(p(1,2))]);
IN_A_Sim(1,:) = [r(1,2) p(1,2)];

[r, p] = corrcoef(Similarity(INs_and_Motifs, 1), Similarity(INs_and_Motifs,3));
disp(['Correlation between IN and Motif: r = ', num2str(r(1,2)), ' and p = ', num2str(p(1,2))]);
IN_Motif_Sim(1,:) = [r(1,2) p(1,2)];

[r, p] = corrcoef(Similarity(INs_and_As, 1), FirstSyllOnsetTime(INs_and_As, 1));
disp(['Correlation between IN similarity and first syll onset time: r = ', num2str(r(1,2)), ' and p = ', num2str(p(1,2))]);
IN_A_OnsetTime(1,:) = [r(1,2) p(1,2)];

disp(['Range of variability in IN similarity is ', num2str(min(Similarity(INs_and_As,1))), ' - ', num2str(max(Similarity(INs_and_As,1)))]);
disp(['Range of variability in first syllable similarity is ', num2str(min(Similarity(INs_and_As,2))), ' - ', num2str(max(Similarity(INs_and_As,2)))]);
disp('Finished analysis of similarity');