function [] = PlotDIRUNDIRSAPFeatures(BatchFileName,Motif)

fid = fopen(BatchFileName,'r');

BatchSAPFeatures = [];

while (~(feof(fid)))
    tline = fgetl(fid);
    NoteFile = [tline,'.not.mat'];
    SAPFeaturesFile = [tline,'.SAPFeatures.mat'];
    
    if ~(exist(SAPFeaturesFile,'file'))
        continue;
    end
    
    if ~(exist(NoteFile,'file'))
        continue;
    end
    
    Notes = load(NoteFile);
    Notes.onsets = Notes.onsets/1000;
    Notes.offsets = Notes.offsets/1000;
    
    SAPFeatures = load(SAPFeaturesFile);
    
    for i = 1:length(Motif),
        MotifNotes = find(Notes.labels == Motif(i));
        for j = 1:length(MotifNotes),
            BatchSAPFeatures(end + 1).Syllable = Motif(i);
            Syllable = (find((SAPFeatures.SAPFeatures.Time >= Notes.onsets(MotifNotes(j))) & (SAPFeatures.SAPFeatures.Time <= Notes.offsets(MotifNotes(j)))));
            BatchSAPFeatures(end).Entropy = mean(SAPFeatures.SAPFeatures.Entropy(Syllable));
            BatchSAPFeatures(end).FM = mean(SAPFeatures.SAPFeatures.FM(Syllable));
            BatchSAPFeatures(end).AM = mean(SAPFeatures.SAPFeatures.AM(Syllable));
            BatchSAPFeatures(end).Amplitude = mean(SAPFeatures.SAPFeatures.Amplitude(Syllable));
            BatchSAPFeatures(end).MeanFreq = mean(SAPFeatures.SAPFeatures.Freq(Syllable));
            BatchSAPFeatures(end).Pitch = mean(SAPFeatures.SAPFeatures.Pitch(Syllable));
            BatchSAPFeatures(end).PitchGoodness = mean(SAPFeatures.SAPFeatures.PitchGoodness(Syllable));
            BatchSAPFeatures(end).Duration = Notes.offsets(MotifNotes(j)) - Notes.onsets(MotifNotes(j));
            if (length(strfind(tline,'_undir')) > 0)
                BatchSAPFeatures(end).DirUndir = 'U';
            else
                BatchSAPFeatures(end).DirUndir = 'D';
            end
        end
    end
end

SAPFeatureNames = fieldnames(BatchSAPFeatures);
SAPFeatureNames = SAPFeatureNames(2:end - 1);
DirSyllables = find([BatchSAPFeatures.DirUndir] == 'D');
UnDirSyllables = find([BatchSAPFeatures.DirUndir] == 'U');

disp(['No of directed syllables = ',num2str(length(DirSyllables))]);
disp(['No of undirected syllables = ',num2str(length(UnDirSyllables))]);

for i = 1:(length(fieldnames(BatchSAPFeatures)) - 3),
    for j = 1:length(BatchSAPFeatures),
        FeatureValues(j) = getfield(BatchSAPFeatures,{j},SAPFeatureNames{i});
    end
    figure;
    set(gcf,'Color','w');
    set(gca,'FontSize',14);
    set(gca,'FontWeight','bold');
    title(SAPFeatureNames{i},'FontSize',16,'FontWeight','bold');
    hold on;
    plot([BatchSAPFeatures(DirSyllables).Duration],FeatureValues(DirSyllables),'b+','MarkerSize',2);
    plot([BatchSAPFeatures(UnDirSyllables).Duration],FeatureValues(UnDirSyllables),'r+','MarkerSize',2);
    xlabel('Duration (sec)');
    ylabel(SAPFeatureNames{i});
    legend('DIR','UNDIR');
end
    
fclose(fid);
