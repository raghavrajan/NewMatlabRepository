function [] = IntroNoteFirstSyllAnalysis(IntroNoteResults, RawDataDir, OutputDir, FileType, OutputFileName, Motif, varargin)

% Using a system for AllINFeatLabels that keeps track of whether the intro
% note was the first, last or a middle intro note. The way I do this is by
% having 3 boolean flags for first, middle and last intro note. For
% instance if there was only one intro note, it would have the flags 1 0 1
% to indicate that is the first, it is also the last.

PresentDir = pwd;

Index = 0;
for i = 1:length(IntroNoteResults.NoofINs),
    Index = Index + 1;
    SongFile = IntroNoteResults.BoutDetails(i).SongFile;
    
    if (strfind(FileType, 'okrank'))
        [TempSong, Fs] = SSAReadOKrankData(RawDataDir, RawDataDir, IntroNoteResults.BoutDetails(i).SongFile, 1);
    else
        if (strfind(FileType, 'wav'))
            cd(RawDataDir);
            [TempSong, Fs] = wavread(SongFile);
            cd(PresentDir);
        else
            if (strfind(FileType, 'obs'))
                channel_string = strcat('obs',num2str(0),'r');
                [TempSong, Fs] = SSASoundIn([RawDataDir, '/'], RawDataDir, SongFile, channel_string);
                % Convert to V - 5V on the data acquisition is 32768
                TempSong = TempSong * 5/32768;
            end
        end
    end
    TempSong = resample(TempSong,44100,Fs);
    Fs = 44100;
    TempSong = TempSong/max(TempSong);
    
    Time = (1:1:length(TempSong))/Fs;
    
    if (IntroNoteResults.NoofINs(i) > 0)
        INs = IntroNoteResults.INs{i};
        INOnset = IntroNoteResults.BoutDetails(i).onsets(INs(end));
        INOffset = IntroNoteResults.BoutDetails(i).offsets(INs(end));
        INWaveform = TempSong(find((Time >= (INOnset - 0.005)) & (Time <= (INOffset + 0.005))));
        cd(OutputDir);
        wavwrite(INWaveform, Fs, 16, [OutputFileName, '.IN.', num2str(Index), '.wav']);
    end
    FirstSyllOnset = IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.MotifStartIndex(i));
    FirstSyllOffset = IntroNoteResults.BoutDetails(i).offsets(IntroNoteResults.MotifStartIndex(i));
    FirstSyllWaveform = TempSong(find((Time >= (FirstSyllOnset - 0.005)) & (Time <= (FirstSyllOffset + 0.005))));
    cd(OutputDir);
    wavwrite(FirstSyllWaveform, Fs, 16, [OutputFileName, '.a.', num2str(Index), '.wav']);
    
    Matches = strfind(IntroNoteResults.BoutDetails(i).labels, Motif);
    if (~isempty(find(Matches == IntroNoteResults.MotifStartIndex(i))))
       MotifOnset = IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.MotifStartIndex(i));
        MotifOffset = IntroNoteResults.BoutDetails(i).offsets(IntroNoteResults.MotifStartIndex(i) + length(Motif) - 1);
        MotifWaveform = TempSong(find((Time >= (MotifOnset - 0.005)) & (Time <= (MotifOffset + 0.005))));
        cd(OutputDir);
        wavwrite(MotifWaveform, Fs, 16, [OutputFileName, '.Motif.', num2str(Index), '.wav']);
    end
end

for i = 1:size(IntroNoteResults.WithinBoutNoofINs,1),
    Index = Index + 1;
    SongFile = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).SongFile;
    if (strfind(FileType, 'okrank'))
        [TempSong, Fs] = SSAReadOKrankData(RawDataDir, RawDataDir, SongFile, 1);
    else
        if (strfind(FileType, 'wav'))
            cd(RawDataDir);
            [TempSong, Fs] = wavread(SongFile);
            cd(PresentDir);
        else
            if (strfind(FileType, 'obs'))
                channel_string = strcat('obs',num2str(0),'r');
                [TempSong, Fs] = SSASoundIn([RawDataDir, '/'], RawDataDir, SongFile, channel_string);
                % Convert to V - 5V on the data acquisition is 32768
                TempSong = TempSong * 5/32768;
            end
        end
    end
    TempSong = resample(TempSong,44100,Fs);
    Fs = 44100;
    TempSong = TempSong/max(TempSong);
    
    Time = (1:1:length(TempSong))/Fs;
    
    if (IntroNoteResults.WithinBoutNoofINs(i,1) > 0)
        INs = IntroNoteResults.WithinBoutINs{i};
        INOnset = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).onsets(INs(end));
        INOffset = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).offsets(INs(end));
        INWaveform = TempSong(find((Time >= (INOnset - 0.005)) & (Time <= (INOffset + 0.005))));
        cd(OutputDir);
        wavwrite(INWaveform, Fs, 16, [OutputFileName, '.IN.', num2str(Index), '.wav']);
    end
    FirstSyllOnset = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).onsets(IntroNoteResults.WithinBoutNoofINs(i, 4));
    FirstSyllOffset = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).offsets(IntroNoteResults.WithinBoutNoofINs(i, 4));
    FirstSyllWaveform = TempSong(find((Time >= (FirstSyllOnset - 0.005)) & (Time <= (FirstSyllOffset + 0.005))));
    cd(OutputDir);
    wavwrite(FirstSyllWaveform, Fs, 16, [OutputFileName, '.a.', num2str(Index), '.wav']);
    
    Matches = strfind(IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).labels, Motif);
    if (~isempty(find(Matches == IntroNoteResults.WithinBoutNoofINs(i, 4))))
        MotifOnset = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).onsets(IntroNoteResults.WithinBoutNoofINs(i, 4));
        MotifOffset = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(i)).offsets(IntroNoteResults.WithinBoutNoofINs(i, 4) + length(Motif) - 1);
        MotifWaveform = TempSong(find((Time >= (MotifOnset - 0.005)) & (Time <= (MotifOffset + 0.005))));
        cd(OutputDir);
        wavwrite(MotifWaveform, Fs, 16, [OutputFileName, '.Motif.', num2str(Index), '.wav']);
    end
end

disp('Finished feature analysis');
