function [] = INAAnalyseSavedNoteAnalysisFiles(FileName)

Temp = load(FileName);
 
figure;
Edges = 0:1:max(Temp.NoofIntroNotes);
Hist = hist(Temp.NoofIntroNotes, Edges);
plot(Edges, Hist/sum(Hist) * 100, 'ks--');
hold on;

Indices = find(Temp.VocalInterval > 0.5);
Hist = hist(Temp.NoofIntroNotes(Indices), Edges);
plot(Edges, Hist/sum(Hist) * 100, 'rs--');

Indices = find(Temp.VocalInterval <= 0.5);
Hist = hist(Temp.NoofIntroNotes(Indices), Edges);
plot(Edges, Hist/sum(Hist) * 100, 'bs--');

Indices = find(Temp.NoofIntroNotes > 0);
NoofINs = Temp.NoofIntroNotes(Indices);
VI = Temp.VocalInterval(Indices);

Hist = hist(NoofINs, Edges);
plot(Edges, Hist/sum(Hist) * 100, 'ks-');

Indices = find(VI > 0.5);
Hist = hist(NoofINs(Indices), Edges);
plot(Edges, Hist/sum(Hist) * 100, 'rs-');

Indices = find(VI <= 0.5);
Hist = hist(NoofINs(Indices), Edges);
plot(Edges, Hist/sum(Hist) * 100, 'bs-');

title(Temp.ANotes.FileName{1});

for i = 1:max(Temp.NoofIntroNotes),
    Freq(i) = length(find(Temp.NoofIntroNotes == i));
end

Indices = find(Freq >= 5);
[MaxVal, MaxInd] = max(Freq);

for i = 1:length(Indices),
    FeatVect{Indices(i)} = [];
    Durs{Indices(i)} = [];
    Ints{Indices(i)} = [];
end

for k = 1:length(Indices),
    for i = 1:length(Temp.MotifIntroNoteIndices),
        if (~isempty(Temp.MotifIntroNoteIndices{i}))
            for j = 1:length(Temp.MotifIntroNoteIndices{i}.NoofINs),
                if (Temp.MotifIntroNoteIndices{i}.NoofINs{j} == Indices(k))
                    FeatVect{Indices(k)} = [FeatVect{Indices(k)}; [Temp.MotifIntroNoteIndices{i}.Entropy{j}' Temp.MotifIntroNoteIndices{i}.Amplitude{j}' Temp.MotifIntroNoteIndices{i}.PG{j}' Temp.MotifIntroNoteIndices{i}.Pitch{j}' Temp.MotifIntroNoteIndices{i}.FM{j}' Temp.MotifIntroNoteIndices{i}.MeanFreq{j}']];
                    Durs{Indices(k)} = [Durs{Indices(k)}; Temp.MotifIntroNoteIndices{i}.Durs{j}'];
                    Ints{Indices(k)} = [Ints{Indices(k)}; Temp.MotifIntroNoteIndices{i}.Intervals{j}'];
                end
            end
        end
    end
end

Mean = mean(Temp.FeatureVect(:,3:8));
STD = std(Temp.FeatureVect(:,3:8));

for k = 1:length(Indices),
    for j = 1:6,
        FeatVect{Indices(k)}(:,j) = (FeatVect{Indices(k)}(:,j) - Mean(j))/STD(j);
    end
end

for k = 1:length(Indices),
    for j = 1:Indices(k),
        TempIndices = j:Indices(k):(size(FeatVect{Indices(k)},1) - Indices(k) + (j));
        MeanFeatVect{Indices(k)}(j,:) = mean(FeatVect{Indices(k)}(TempIndices,:));
    end
end

Cols = ['rgbcmk'];
figure;
hold on;
for k = 1:length(Indices),
        plot(mean(Durs{Indices(k)}), mean(Ints{Indices(k)}), [Cols(k), 's-']);
end
title(Temp.ANotes.FileName{1});

figure;
hold on;
for k = 1:length(Indices),
    MeanDistVect = [MeanFeatVect{MaxInd}(end,:)];
    D = squareform(pdist([MeanDistVect; MeanFeatVect{Indices(k)}]));
    plot([-(Indices(k) - 1):1:0], D(1,(end-Indices(k) + 1):end), [Cols(k), 's-']);
end
title(Temp.ANotes.FileName{1});

figure;
hold on;
Cols = ['rgbcmk'];
for k = 1:length(Indices),
    MeanDistVect = [MeanFeatVect{MaxInd}(1,:)];
    D = squareform(pdist([MeanDistVect; MeanFeatVect{Indices(k)}]));
    plot([-(Indices(k) - 1):1:0], D(1,(end-Indices(k) + 1):end), [Cols(k), 's-']);
end
title(Temp.ANotes.FileName{1});