function [] = MakeManySpectrogramsCoverFigure(DataDir, FileList, FileType, SpectColourMap, NumBouts, Motif)

OutputDir = '/home/raghav/StudentRelated/Divya/FirstManuscript';
PresentDir = pwd;
cd(DataDir);

Fid = fopen(FileList, 'r');
SongFiles = textscan(Fid, '%s', 'DeLimiter', '\n');
SongFiles = SongFiles{1};
fclose(Fid);


% Now load up note files
for i = 1:length(SongFiles),
    NoteInfo{i} = load(fullfile('ASSLNoteFiles', [SongFiles{i}, '.not.mat']));
end

% Now randomly find files with valid song bouts and then plot spectrograms
% of them with a random interval at the beginning

figure;
set(gcf, 'Color', 'k');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 23 24]);
set(gcf, 'ReSize', 'off');
p = panel();
p.pack({1.5/7 1/7 4.5/7});
p(1).pack(2*NumBouts/4, 1);
p(3).pack(3*(2*NumBouts/4), 1);
p.de.margin = 1;

InterBoutInterval = 1950; % in ms
Index = 1;
i = 0;
while ((Index <= (2*NumBouts + 1)) && (i < length(SongFiles)))
    i = i + 1;
    if (isempty(strfind(NoteInfo{i}.labels, Motif)))
        continue;
    end
    [RawData, Fs] = GetData(DataDir, SongFiles{i}, FileType, 0);
    
    % Now find bouts and pick out a song bout and plot spectrogram of it
    Intervals = [NoteInfo{i}.onsets(:); length(RawData)*1000/Fs] - [0; NoteInfo{i}.offsets(:)];
    Bouts = find(Intervals >= InterBoutInterval);
    
    if (length(Bouts) < 2)
        continue;
    end
    
    for j = 1:length(Bouts)-1,
        BoutOnsetTime = NoteInfo{i}.onsets(Bouts(j))/1000 - 0.2;
        BoutOffsetTime = NoteInfo{i}.offsets(Bouts(j+1)-1)/1000;
        
        if (isempty(strfind(NoteInfo{i}.labels(Bouts(j):(Bouts(j+1)-1)), Motif)))
            continue;
        else
            Motifs = (Bouts(j) - 1) + strfind(NoteInfo{i}.labels(Bouts(j):(Bouts(j+1)-1)), Motif);
            MotifOnsetTime = NoteInfo{i}.onsets(Motifs(1))/1000;
            MotifOffsetTime = 0.075 + NoteInfo{i}.offsets(Motifs(1)+length(Motif)-1)/1000;
            MotifLen(Index) = MotifOffsetTime - MotifOnsetTime;
        end
        
        if (Index <= (2*NumBouts/4))
            p(1,Index,1).select();
        else
            if (Index == ((2*NumBouts/4)+1))
                p(2).select();
            else
                p(3,Index-((2*NumBouts/4)+1),1).select();
            end
        end
        
        PlotSpectrogramInAxis_SongVar_ColorMap(RawData(round((BoutOnsetTime) * Fs):round((MotifOffsetTime)*Fs)), linspace(BoutOnsetTime - MotifOnsetTime, MotifOffsetTime - MotifOnsetTime, length(RawData(round((BoutOnsetTime) * Fs):round((MotifOffsetTime)*Fs)))), Fs, gca, SpectColourMap);
        set(gca, 'Visible', 'off');
        Index = Index + 1;
        if (Index > (2*NumBouts+1))
            break;
        end
    end
end
p(2).select();
axis([-1.5 max(MotifLen)+0.05 300 8000]);

for Index = 1:2*NumBouts/4,
    p(1,Index, 1).select();
    axis([-1.5 max(MotifLen)+0.05 300 8000]);
end
for Index = 1:3*2*NumBouts/4,
    p(3,Index, 1).select();
    axis([-1.5 max(MotifLen)+0.05 300 8000]);
end
%map = colormap(SpectColourMap);
% colormap(1 - map);
p.marginbottom = 1;
p.marginleft = 1;
p.marginright = 1;
p.margintop = 1;
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'InvertHardcopy', 'off');
print(fullfile(OutputDir, 'CoverMultipleSpectrograms.png'), '-r300', '-dpng');
