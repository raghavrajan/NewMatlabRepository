function [] = PlotSyllableTrajectories(DirectoryName, FileList, NoteFileDir, FileType)

PresentDir = pwd;
cd(DirectoryName);

% First get all the files
Fid = fopen(FileList, 'r');
Files = textscan(Fid, '%s', 'DeLimiter', '\n');
Files = Files{1};
fclose(Fid);

% Now load up the note files
for i = 1:length(Files),
    NoteInfo{i} = load(fullfile(DirectoryName, NoteFileDir, [Files{i}, '.not.mat']));
end

% Now put together spectrograms of multiple continuous files together and
% then do a pca on this 
CombinedSpect = [];
for i = 1:min(10,length(Files)),
    [RawData, Fs] = GetData(DirectoryName, Files{i}, FileType, 0);
    [PowSpect, freq, spect, time_song] = CalculateMultiTaperSpectrogram(RawData, Fs, 8, 4, 1.5);
    Freq = find((freq >= 860) & (freq <= 8000));
    CombinedSpect = [CombinedSpect PowSpect(Freq,:)];
end

% Now do pca with frequencies as the variables
[Coeff, Score, Latent] = pca(CombinedSpect');

disp('Finished');