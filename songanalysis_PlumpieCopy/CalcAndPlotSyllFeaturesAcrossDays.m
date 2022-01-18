function [] = CalcAndPlotSyllFeaturesAcrossDays(DirectoryList, FileList, FileType)

PresentDir = pwd;
Index = 0;
for i = 1:length(DirectoryList),
    fprintf(' >> %d ', i);
    Fid = fopen(fullfile(DirectoryList{i}, FileList{i}), 'r');
    Files = textscan(Fid, '%s', 'DeLimiter', '\n');
    Files = Files{1};
    fclose(Fid);
    for j = 1:length(Files),
        Index = Index + 1;
        AllDirectories{Index} = DirectoryList{i};
        AllFiles{Index} = Files{j};
        [RawData, Fs] = GetData(DirectoryList{i}, Files{j}, FileType, 1);
        AllNotes{Index} = load(fullfile(DirectoryList{i}, 'ASSLNoteFiles', [Files{j}, '.not.mat']));
        [Feats, RawFeats, FeatsFs] = ASSLCalculateSAPFeatsWithOnsets_SplitSylls(RawData, (1:1:length(RawData))/Fs, Fs, AllNotes{Index}.onsets/1000, AllNotes{Index}.offsets/1000);
        AllFeats{Index} = Feats;
    end
end
fprintf('\n');
disp('Finished');