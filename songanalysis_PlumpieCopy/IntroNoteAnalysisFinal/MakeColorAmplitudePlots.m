function [] = MakeColorAmplitudePlots(FileList, Motif, NoteFileDir)

Fid = fopen(FileList, 'r');
Temp = textscan(Fid, '%s', 'delimiter', '\n');
fclose(Fid);

SongFiles = Temp{1};
Resolution = 0.001; % seconds
TimeAxis = -2.5:Resolution:5;

Amplitude = ones(length(SongFiles), length(TimeAxis)) * -0.5;

NoOfIterations = min(25, length(SongFiles));

for i = 1:NoOfIterations,
    SlashIndex = find((SongFiles{i} == '/') | (SongFiles{i} == '\'));
    if (~isempty(SlashIndex))
        SongFiles{i} = SongFiles{i}(SlashIndex(end)+1:end);
    end
    Notes = load([NoteFileDir, '/', SongFiles{i}, '.not.mat']);
    Notes.onsets = Notes.onsets/1000;
    Notes.offsets = Notes.offsets/1000;
    FirstMotifSyll = find(Notes.labels == Motif(1));
    if (~isempty(FirstMotifSyll))
        FirstMotifSyll = FirstMotifSyll(1);
        FirstMotifSyllOnsetTime = Notes.onsets(FirstMotifSyll);
        Notes.onsets = Notes.onsets - FirstMotifSyllOnsetTime;
        Notes.offsets = Notes.offsets - FirstMotifSyllOnsetTime;
        Sylls = find((Notes.onsets >= TimeAxis(1)) & (Notes.offsets <= TimeAxis(end)));
        if (~isempty(Sylls))
            Amplitude(i, 1:round((Notes.onsets(Sylls(1)) - TimeAxis(1))/Resolution)) = 0;
            for j = 1:length(Sylls),
                if (~isempty(find(Motif == Notes.labels(Sylls(j)))))
                    Amplitude(i,round((Notes.onsets(Sylls(j)) - TimeAxis(1))/Resolution):round((Notes.offsets(Sylls(j)) - TimeAxis(1))/Resolution)) = 1;
                else
                    Amplitude(i,round((Notes.onsets(Sylls(j)) - TimeAxis(1))/Resolution):round((Notes.offsets(Sylls(j)) - TimeAxis(1))/Resolution)) = 0.5;
                end
                
                if (j ~= length(Sylls))
                    Amplitude(i,round((Notes.offsets(Sylls(j)) - TimeAxis(1))/Resolution):round((Notes.onsets(Sylls(j+1)) - TimeAxis(1))/Resolution)) = 0;
                else
                    Amplitude(i,round((Notes.offsets(Sylls(j)) - TimeAxis(1))/Resolution):end) = 0;
                end
            end
        end
    end
end

figure;
imagesc(TimeAxis, 1:1:length(SongFiles), Amplitude);
axis([-1.5 0.5 0.5 (NoOfIterations + 0.5)]);

disp('Finished');
       