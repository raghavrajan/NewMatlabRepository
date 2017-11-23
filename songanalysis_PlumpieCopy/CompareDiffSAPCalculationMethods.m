function [] = CompareDiffSAPCalculationMethods(DataDir, SongFile, FileType, NoteFileDir)

[RawData, Fs] = GetData(DataDir, SongFile, FileType, 1);

NoteInfo = load(fullfile(DataDir, NoteFileDir, [SongFile, '.not.mat']));

Time = (1:1:length(RawData))/Fs;

tic;
[Feats, RawFeats, FeatsFs] = ASSLCalculateSAPFeatsWithOnsets(RawData, Time, Fs, NoteInfo.onsets/1000, NoteInfo.offsets/1000);
toc

tic;
[Feats2, RawFeats2, FeatsFs2] = ASSLCalculateSAPFeatsWithOnsets_SplitSylls(RawData, Time, Fs, NoteInfo.onsets/1000, NoteInfo.offsets/1000);
toc

FeatFields = fieldnames(Feats);

for i = 1:length(FeatFields),
    [min(eval(['Feats.', FeatFields{i}]) - eval(['Feats2.', FeatFields{i}])) max(eval(['Feats.', FeatFields{i}]) - eval(['Feats2.', FeatFields{i}]))]
end

disp('Finished');


