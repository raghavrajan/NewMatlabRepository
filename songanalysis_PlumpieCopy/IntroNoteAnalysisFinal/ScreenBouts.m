function [BoutDetails, AllBoutDetails] = ScreenBouts(RawDataDir, RecFileDir, FileType, NoteFileDir, SongFileList, InterBoutInterval, MotifSylls, ContinuousOrDiscrete, varargin)

PresentDir = pwd;
if (strfind(ContinuousOrDiscrete, 'continuous'))
    ContinuousFileTime = varargin{1};
    if (length(varargin) > 1)
        ContinuousFiles = varargin{2};
    end
end

TotalBouts = 0;
ShortStarts = 0;
ShortEnds = 0;
RealBouts = 0;
AllBouts = 0;

Fid = fopen(SongFileList, 'r');
TempNoteFiles = textscan(Fid, '%s', 'DeLimiter', '\n');
NoteFiles = TempNoteFiles{1};
fclose(Fid);

if (strfind(ContinuousOrDiscrete, 'continuous'))
    if (~exist('ContinuousFiles', 'var'))
        for i = 1:length(NoteFiles),    
            FileTimes(i) = str2double(NoteFiles{i}(end-5:end-4))*3600 + str2double(NoteFiles{i}(end-3:end-2))*60 + str2double(NoteFiles{i}(end-1:end));
        end
        Diff_FileTimes = diff(FileTimes);
        BreakIndices = find(Diff_FileTimes > ContinuousFileTime);
        if (~isempty(BreakIndices))
            ContinuousFiles{1} = [1:BreakIndices(1)];
    
            for i = 1:length(BreakIndices)-1,
                ContinuousFiles{i+1} = [BreakIndices(i)+1:BreakIndices(i+1)];
            end
            ContinuousFiles{end+1} = [BreakIndices(end)+1:length(NoteFiles)];
        else
            ContinuousFiles{1} = [1:length(FileTimes)];
        end
    end
end

if (strfind(ContinuousOrDiscrete, 'continuous'))
    for FileNo = 1:length(ContinuousFiles),
        Song = [];
        ActualNotes.onsets = [];
        ActualNotes.offsets = [];
        ActualNotes.labels = [];
        ActualNotes.Fs = [];
        TimeOffset = 0;
        for i = ContinuousFiles{FileNo},
            SlashIndex = find((NoteFiles{i} == '/') | (NoteFiles{i} == '\'));
            if (~isempty(SlashIndex))
                SongFile = NoteFiles{i}(SlashIndex(end)+1:end);
            else
                SongFile = NoteFiles{i};
            end

            disp(SongFile);

            if (strfind(FileType, 'okrank'))
                [TempSong, Fs] = SSAReadOKrankData(RawDataDir, RecFileDir, SongFile, 1);
                Song = [Song; TempSong];
            else
                if (strfind(FileType, 'wav'))
                    cd(RawDataDir);
                    [TempSong, Fs] = wavread(SongFile);
                    Song = [Song; TempSong];
                    cd(PresentDir);
                else
                    if (strfind(FileType, 'obs'))
                        channel_string = strcat('obs',num2str(0),'r');
                        [TempSong, Fs] = SSASoundIn([RawDataDir, '/'], RecFileDir, SongFile, channel_string);
                        % Convert to V - 5V on the data acquisition is 32768
                        TempSong = TempSong * 5/32768;
                        Song = [Song; TempSong];
                    end
                end
            end

            cd(NoteFileDir);
            Notes = load([NoteFiles{i}, '.not.mat']);
            Notes.onsets = Notes.onsets(:);
            Notes.offsets = Notes.offsets(:);
            
            Notes.onsets(find(Notes.labels == '0')) = [];
            Notes.offsets(find(Notes.labels == '0')) = [];
            Notes.labels(find(Notes.labels == '0')) = [];

            if (isempty(Notes.onsets))
                TimeOffset = TimeOffset + length(TempSong)/Fs;            
                continue;
            else
                ActualNotes.onsets = [ActualNotes.onsets; (Notes.onsets + TimeOffset*1000)];
                ActualNotes.offsets = [ActualNotes.offsets; (Notes.offsets + TimeOffset*1000)];
                ActualNotes.labels = [ActualNotes.labels Notes.labels];
                ActualNotes.Fs = [ActualNotes.Fs; Notes.Fs];
                TimeOffset = TimeOffset + length(TempSong)/Fs;            
            end
            cd(PresentDir);
        end
        clear Notes;
        if (isempty(ActualNotes.onsets))
            continue;
        end
        
        Time = (1:1:length(Song))/Fs;
        
        RepeatedSylls = regexp(ActualNotes.labels, '[A-Z]');
        if (~isempty(RepeatedSylls))
            Counter = 1;
            SyllsToBeChanged = [];
            while (Counter < (length(RepeatedSylls)))
                if (RepeatedSylls(Counter) == (RepeatedSylls(Counter+1) - 1))
                    if (ActualNotes.labels(RepeatedSylls(Counter)) == ActualNotes.labels(RepeatedSylls(Counter+1)))
                        SyllsToBeChanged = [SyllsToBeChanged; RepeatedSylls(Counter)];
                        Counter = Counter + 2;
                    else
                        Counter = Counter + 1;
                    end
                else
                    Counter = Counter + 1;
                end
            end
            ActualNotes.labels(SyllsToBeChanged) = lower(ActualNotes.labels(SyllsToBeChanged));
            ActualNotes.labels(SyllsToBeChanged+1) = [];
            ActualNotes.onsets(SyllsToBeChanged+1) = [];
            ActualNotes.offsets(SyllsToBeChanged) = [];
        end

        Notes = ActualNotes;
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
                AllBouts = AllBouts + 1;
                AllBoutDetails(AllBouts).SongFile = NoteFiles(ContinuousFiles{FileNo});
                AllBoutDetails(AllBouts).Fs = Fs;
                AllBoutDetails(AllBouts).BoutIndices = Bouts(j,:);
                AllBoutDetails(AllBouts).BoutOnset = Notes.onsets(Bouts(j,1))/1000 - InterBoutInterval;
                AllBoutDetails(AllBouts).BoutOffset = Notes.offsets(Bouts(j,2))/1000 + InterBoutInterval;
                AllBoutDetails(AllBouts).BoutLength = BoutLengths(j);
                AllBoutDetails(AllBouts).onsets = Notes.onsets(Bouts(j,1):1:Bouts(j,2))/1000;
                AllBoutDetails(AllBouts).offsets = Notes.offsets(Bouts(j,1):1:Bouts(j,2))/1000;
                AllBoutDetails(AllBouts).labels = Notes.labels(Bouts(j,1):1:Bouts(j,2));

                %for Note = 1:length(AllBoutDetails(AllBouts).onsets),
                %    SyllStartIndex = find(Time <= (AllBoutDetails(AllBouts).onsets(Note)), 1, 'last') - round(0.05*Fs);
                %    SyllEndIndex = find(Time <= (AllBoutDetails(AllBouts).offsets(Note)), 1, 'last') + round(0.05*Fs);
                %    AllBoutDetails(AllBouts).FFTAmplitudeWaveform{Note} = smooth(SyllStartIndex:SyllEndIndex);
                %    %[m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude ,m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ]=deriv(Song(SyllStartIndex:SyllEndIndex), Fs);
                    %m_amplitudeTime = linspace(Time(SyllStartIndex), Time(SyllEndIndex), length(m_amplitude));
                    %AllBoutDetails(AllBouts).SAPAmplitudeWaveform{Note} = spline(m_amplitudeTime, m_amplitude, Time(SyllStartIndex:SyllEndIndex));
                    %AllBoutDetails(AllBouts).AmplitudeWaveformStartEndTime(Note,:) = [-0.05 0.05];
                %end
                StartIndex = round(AllBoutDetails(AllBouts).BoutOnset * Fs);
                EndIndex = round(AllBoutDetails(AllBouts).BoutOffset * Fs);

                AllBoutDetails(AllBouts).Feats = CalculateSAPFeatsWithOnsets(Song(StartIndex:EndIndex), Time(StartIndex:EndIndex), Fs, AllBoutDetails(AllBouts).onsets, AllBoutDetails(AllBouts).offsets);
                
                SyllFlag = ones(size(MotifSylls));
                for k = 1:length(MotifSylls),
                    if (isempty(find(Notes.labels(Bouts(j,1):1:Bouts(j,2)) == MotifSylls(k))))
                        SyllFlag(k) = 0;
                    end
                end

                if (~isempty(find(SyllFlag == 1)))
                    AllBoutDetails(AllBouts).Motif = 1;
                    RealBouts = RealBouts + 1;
                    BoutDetails(RealBouts).SongFile = NoteFiles(ContinuousFiles{FileNo});
                    BoutDetails(RealBouts).Fs = Fs;
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
                else
                    AllBoutDetails(AllBouts).Motif = 0;
                end
            end
        end
    end
else
    for i = 1:length(NoteFiles),
        cd(NoteFileDir);
        NoteFiles{i} = [NoteFiles{i}, '.not.mat'];
        Notes = load(NoteFiles{i});
        Notes.onsets = Notes.onsets(:);
        Notes.offsets = Notes.offsets(:);
        
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
                AllBouts = AllBouts + 1;
                AllBoutDetails(AllBouts).SongFile = SongFile;
                AllBoutDetails(AllBouts).Fs = Fs;
                AllBoutDetails(AllBouts).BoutIndices = Bouts(j,:);
                AllBoutDetails(AllBouts).BoutOnset = Notes.onsets(Bouts(j,1))/1000 - InterBoutInterval;
                AllBoutDetails(AllBouts).BoutOffset = Notes.offsets(Bouts(j,2))/1000 + InterBoutInterval;
                AllBoutDetails(AllBouts).BoutLength = BoutLengths(j);
                AllBoutDetails(AllBouts).onsets = Notes.onsets(Bouts(j,1):1:Bouts(j,2))/1000;
                AllBoutDetails(AllBouts).offsets = Notes.offsets(Bouts(j,1):1:Bouts(j,2))/1000;
                AllBoutDetails(AllBouts).labels = Notes.labels(Bouts(j,1):1:Bouts(j,2));

                StartIndex = round(AllBoutDetails(AllBouts).BoutOnset * Fs);
                EndIndex = round(AllBoutDetails(AllBouts).BoutOffset * Fs);

                AllBoutDetails(AllBouts).Feats = CalculateSAPFeatsWithOnsets(Song(StartIndex:EndIndex), Time(StartIndex:EndIndex), Fs, AllBoutDetails(AllBouts).onsets, AllBoutDetails(AllBouts).offsets);
                
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
end

disp(['Total Bouts in the ', num2str(length(NoteFiles)), ' note files: ', num2str(TotalBouts)]);
disp(['Number of bouts with inter bout intervals greater than ', num2str(InterBoutInterval), 's: ', num2str(RealBouts)]);
