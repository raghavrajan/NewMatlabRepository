function [] = SeparateSyllables(BatchFileName,Motif)

fid = fopen(BatchFileName,'r');
DirectoryName = pwd;
DirectoryName(end + 1) = '/';

Silence = wavread('Silence.wav');

while (~(feof(fid)))
    tline = fgetl(fid);
    NoteFile = [tline,'.not.mat'];
    FiltFile = [tline,'.filt'];
    
    if ~(exist(NoteFile,'file'))
        continue;
    end
    
    Notes = load(NoteFile);
    
    [FiltSong, Fs] = soundin_copy(DirectoryName,FiltFile,'filt');
    FiltSong = FiltSong/32768;
    
    TimeSong = 0:1/Fs:(length(FiltSong) - 1)/Fs;
    TimeSong = TimeSong * 1000; % convert time to ms
    
    for i = 1:length(Motif),
        MotifNotes = find(Notes.labels == Motif(i));
        for j = 1:length(MotifNotes),
            SyllableIndices(1) = find(TimeSong <= (Notes.onsets(MotifNotes(j))),1,'last');
            SyllableIndices(2) = find(TimeSong <= (Notes.offsets(MotifNotes(j))),1,'last');
            FiltSyllable = FiltSong(SyllableIndices(1):SyllableIndices(2));
            FiltSyllable = [Silence; FiltSyllable; Silence];
            TimeSyllable = 0:1/Fs:(length(FiltSyllable) - 1)/Fs;
            plot(TimeSyllable,FiltSyllable);
            cd SeparateSyllables;
            OutputFileName = [tline,'_',Motif(i),'_',num2str(j),'.wav'];
            wavwrite(FiltSyllable,44100,16,OutputFileName);
            cd ../;
        end
    end
end

fclose(fid);
