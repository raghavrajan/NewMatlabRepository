function [] = PlotSpectrograms(DirectoryName, FileName, FileExt, FileNos, FileType, Motif, TitleString)

PresentDirectory = pwd;
cd(DirectoryName);

FigureIndex = 1;
SongSpectrogramFigure(FigureIndex) = figure;
hold on;
SongIndex = 0;

if (length(FileNos.Directed) > 0)
    for i = 1:length(FileNos.Directed),
        if (strfind(FileType, 'okrank'))
            if (FileNos.Directed(i) > 100000)
                DataFileName = [FileName,num2str(FileNos.Directed(i))];
            else
                DataFileName = [FileName,'0',num2str(FileNos.Directed(i))];
            end
            NoteFileName = [DataFileName,'.not.mat'];
            [Data, Fs] = ReadOKrankData(DirectoryName, DataFileName, 0);
        end
        
        Notes = load(NoteFileName);
        Motifs = strfind(Notes.labels, Motif);
        if (length(Motifs) > 0)
            for i = 1:length(Motifs),
                StartIndex = round(Notes.onsets(Motifs(i)) * Fs/1000);
                EndIndex = round(Notes.offsets(Motifs(i) + length(Motif) - 1) * Fs/1000);
                RawSong = Data(StartIndex:EndIndex);
                
                Time = 0:1/Fs:(length(RawSong) - 1)/Fs;
                
                % Now using an 8 pole butterworth bandpass filter as default.
                [b,a]=butter(8,[300*2/Fs, 10000*2/Fs]);

                FiltSong=filtfilt(b, a, RawSong);
  
                if (length(RawSong) ~= length(FiltSong))
                    disp(['warning! bandpass: input and output file lengths do not match!']);
                end


                nfft=round(Fs*8/1000);
                nfft = 2^nextpow2(nfft);
                spect_win = hanning(nfft);
                noverlap = round(0.9*length(spect_win)); %number of overlapping points       
                %now calculate spectrogram
                [spect, freq, time] = spectrogram(FiltSong,spect_win,noverlap,nfft,Fs,'yaxis');
                freq = freq + (SongIndex * 12000);
                idx_spect=scale_spect(spect);  %calculate index array for spectrogram
                f_min = freq(1);
                f_max = freq(length(freq));
                freq_spect = [f_min, f_max];
                t_min = time(1);
                t_max = time(end);

                %adjust time axis for spectrogram offset (1/2 window duration in ms)
                t_min = t_min + 0.5*nfft/Fs;  
                t_max = t_max + 0.5*nfft/Fs;  

                time_spect = [t_min, t_max];                
                
                figure(SongSpectrogramFigure(FigureIndex));
                hold on;

                cm = disp_idx_spect(idx_spect, time_spect, freq_spect, -55, 0, 1.2, 'gray', 'classic');
                
                SongIndex = SongIndex + 1;
                
                if (SongIndex > 10)
                    FigureIndex = FigureIndex + 1;
                    SongSpectrogramFigure(FigureIndex) = figure;
                    hold on;
                end
            end
        end
    end
end
   
if (length(FileNos.Undirected) > 0)
    for i = 1:length(FileNos.Undirected),
        if (strfind(FileType, 'okrank'))
            if (FileNos.Undirected(i) > 100000)
                DataFileName = [FileName,num2str(FileNos.Undirected(i))];
            else
                DataFileName = [FileName,'0',num2str(FileNos.Undirected(i))];
            end
            NoteFileName = [DataFileName,'.not.mat'];
            [Data, Fs] = ReadOKrankData(DirectoryName, DataFileName, 0);
        end
        
        Notes = load(NoteFileName);
        Motifs = strfind(Notes.labels, Motif);
        if (length(Motifs) > 0)
            for i = 1:length(Motifs),
                StartIndex = round(Notes.onsets(Motifs(i)) * Fs/1000);
                EndIndex = round(Notes.offsets(Motifs(i) + length(Motif) - 1) * Fs/1000);
                RawSong = Data(StartIndex:EndIndex);
                
                Time = 0:1/Fs:(length(RawSong) - 1)/Fs;
                
                % Now using an 8 pole butterworth bandpass filter as default.
                [b,a]=butter(8,[300*2/Fs, 10000*2/Fs]);

                FiltSong=filtfilt(b, a, RawSong);
  
                if (length(RawSong) ~= length(FiltSong))
                    disp(['warning! bandpass: input and output file lengths do not match!']);
                end


                nfft=round(Fs*8/1000);
                nfft = 2^nextpow2(nfft);
                spect_win = hanning(nfft);
                noverlap = round(0.9*length(spect_win)); %number of overlapping points       
                %now calculate spectrogram
                [spect, freq, time] = spectrogram(FiltSong,spect_win,noverlap,nfft,Fs,'yaxis');
                freq = freq + (SongIndex * 12000);
                idx_spect=scale_spect(spect);  %calculate index array for spectrogram
                f_min = freq(1);
                f_max = freq(length(freq));
                freq_spect = [f_min, f_max];
                t_min = time(1);
                t_max = time(end);

                %adjust time axis for spectrogram offset (1/2 window duration in ms)
                t_min = t_min + 0.5*nfft/Fs;  
                t_max = t_max + 0.5*nfft/Fs;  

                time_spect = [t_min, t_max];                
                
                figure(SongSpectrogramFigure(FigureIndex));
                hold on;

                if (~exist('cm','var'))
                    cm = disp_idx_spect(idx_spect, time_spect, freq_spect, -55, 0, 1.0, 'gray', 'classic');
                end
                image(time_spect, freq_spect, idx_spect)
                colormap(cm);
                
                SongIndex = SongIndex + 1;
                
                if (SongIndex > 10)
                    FigureIndex = FigureIndex + 1;
                    SongSpectrogramFigure(FigureIndex) = figure;
                    SongIndex = 0;
                    hold on;
                end
            end
        end
    end
end

for i = 1:FigureIndex,
    figure(SongSpectrogramFigure(i));
    set(gcf,'Color','w');
    title(TitleString);
    axis tight;
end
