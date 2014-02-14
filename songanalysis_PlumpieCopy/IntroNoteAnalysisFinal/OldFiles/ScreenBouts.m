function [BoutDetails] = ScreenBouts(RawDataDir, RecFileDir, FileType, NoteFileDir, SongFileList, InterBoutInterval, MotifSylls)

PresentDir = pwd;

TotalBouts = 0;
ShortStarts = 0;
ShortEnds = 0;
RealBouts = 0;

Fid = fopen(SongFileList, 'r');
TempNoteFiles = textscan(Fid, '%s', 'DeLimiter', '\n');
NoteFiles = TempNoteFiles{1};
fclose(Fid);

for i = 1:length(NoteFiles),
    cd(NoteFileDir);
    NoteFiles{i} = [NoteFiles{i}, '.not.mat'];
    Notes = load(NoteFiles{i});
    cd(PresentDir);
    
    Notes.onsets(find(Notes.labels == '0')) = [];
    Notes.offsets(find(Notes.labels == '0')) = [];
    Notes.labels(find(Notes.labels == '0')) = [];
    
    if (isempty(Notes.onsets))
        continue;
    end
    
    SlashIndex = find((NoteFiles{i} == '/') | (NoteFiles{i} == '\'));
    if (~isempty(SlashIndex))
        SongFile = NoteFiles{i}(SlashIndex(end)+1:end-8);
    else
        SongFile = NoteFiles{i}(1:end-8);
    end
    
    disp(SongFile);
    
    if (strfind(FileType, 'okrank'))
        [Song, Fs] = SSAReadOKrankData(RawDataDir, RecFileDir, SongFile, 1);
    else
        if (strfind(FileType, 'wav'))
            cd(RawDataDir);
            [Song, Fs] = wavread(SongFile);
            cd(PresentDir);
        else
            if (strfind(FileType, 'obs'))
                channel_string = strcat('obs',num2str(0),'r');
                [Song, Fs] = SSASoundIn([RawDataDir, '/'], RecFileDir, SongFile, channel_string);

                % Convert to V - 5V on the data acquisition is 32768
                Song = Song * 5/32768;
            end
        end
    end
    
    Time = (1:1:length(Song))/Fs;
    cd(PresentDir);
    
    FileLength = length(Song)/Fs;

    Intervals = Notes.onsets(2:end)/1000 - Notes.offsets(1:end-1)/1000;
    TempBouts = find(Intervals >= InterBoutInterval);
    Bouts = [];
    BoutLengths = [];
    
    for j = 0:length(TempBouts),
        if (j ~= length(TempBouts))
            if (j == 0)
                Bouts = [Bouts; [1 TempBouts(j+1)]];
                BoutLengths = [BoutLengths; (Notes.offsets(TempBouts(j+1)) - Notes.onsets(1))/1000];
            else
                Bouts = [Bouts; [(TempBouts(j) + 1) TempBouts(j+1)]];
                BoutLengths = [BoutLengths; (Notes.offsets(TempBouts(j+1)) - Notes.onsets(TempBouts(j) + 1))/1000];
            end
        else
            if (j == 0)
                Bouts = [Bouts; [1 length(Notes.onsets)]];
                BoutLengths = [BoutLengths; (Notes.offsets(end) - Notes.onsets(1))/1000];
            else
                Bouts = [Bouts; [(TempBouts(j) + 1) length(Notes.onsets)]];
                BoutLengths = [BoutLengths; (Notes.offsets(end) - Notes.onsets(TempBouts(j) + 1))/1000];
            end
        end
    end
    
   
    for j = 1:size(Bouts,1),
        Flag = 1;
        TotalBouts = TotalBouts + 1;
        if (Notes.onsets(Bouts(j,1))/1000 <= InterBoutInterval)
            ShortStarts = ShortStarts + 1;
            Flag = 0;
        end
        
        if (Notes.offsets(Bouts(j,2))/1000 > (FileLength - InterBoutInterval))
            ShortEnds = ShortEnds + 1;
            Flag = 0;
        end
        
        if (Flag == 1)
            SyllFlag = ones(size(MotifSylls));
            for k = 1:length(MotifSylls),
                if (isempty(find(Notes.labels(Bouts(j,1):1:Bouts(j,2)) == MotifSylls(k))))
                    SyllFlag(k) = 0;
                end
            end
            
            if (~isempty(find(SyllFlag == 1)))
                RealBouts = RealBouts + 1;
                BoutDetails(RealBouts).SongFile = SongFile;
                BoutDetails(RealBouts).BoutIndices = Bouts(j,:);
                BoutDetails(RealBouts).BoutOnset = Notes.onsets(Bouts(j,1))/1000 - InterBoutInterval;
                BoutDetails(RealBouts).BoutOffset = Notes.offsets(Bouts(j,2))/1000 + InterBoutInterval;
                BoutDetails(RealBouts).BoutLength = BoutLengths(j);
                BoutDetails(RealBouts).onsets = Notes.onsets(Bouts(j,1):1:Bouts(j,2))/1000;
                BoutDetails(RealBouts).offsets = Notes.offsets(Bouts(j,1):1:Bouts(j,2))/1000;
                BoutDetails(RealBouts).labels = Notes.labels(Bouts(j,1):1:Bouts(j,2));
                
                StartIndex = round(BoutDetails(RealBouts).BoutOnset * Fs);
                EndIndex = round(BoutDetails(RealBouts).BoutOffset * Fs);
                
                BoutDetails(RealBouts).Feats = CalculateSAPFeatsWithOnsets(Song(StartIndex:EndIndex), Time(StartIndex:EndIndex), Fs, BoutDetails(RealBouts).onsets, BoutDetails(RealBouts).offsets);
            end
        end
    end
end

disp(['Total Bouts in the ', num2str(length(NoteFiles)), ' note files: ', num2str(TotalBouts)]);
disp(['Number of bouts with inter bout intervals greater than ', num2str(InterBoutInterval), 's: ', num2str(RealBouts)]);
