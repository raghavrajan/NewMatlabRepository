function [Bouts] = JSongAnalysisBatch(DirectoryName, BirdName, FileExt, FileType, RMS)

if (DirectoryName(end) ~= '/')
    DirectoryName(end + 1) = '/';
end
    
cd(DirectoryName);

DataFiles = dir([BirdName,'*',FileExt]);

BoutIndex = 1;
Bouts = [];

if (strfind(FileType, 'obs'))
    for i = 1:1000,
        try
            [RawSong, Fs] = soundin_copy(DirectoryName, DataFiles(i).name, FileType);
            RawSong = RawSong/32768 * 1;
            FiltSong = bandpass_fft_filter(RawSong, 500, 100, Fs);
%           FiltSong = bandpass(RawSong, Fs, 300, 8000, 'butter');
            SquaredSong = FiltSong.^2;
            Len = round(Fs*2/1000);
            h = ones(1, Len)/Len;
            Smooth = conv(h, SquaredSong);
            Offset = round((length(Smooth) - length(FiltSong))/2);
            Smooth = Smooth((1 + Offset):(length(FiltSong) + Offset));

            NoteTimes = Smooth > 10*RMS;
            
            Onsets = [];
            Offsets = [];
            Transitions = find(diff(NoteTimes));
            
            if (length(NoteTimes > 1))
            
                if (size(NoteTimes,1) == 1)
                    NoteTimes = NoteTimes';
                end

                if (NoteTimes(1) == 1)
                    if (NoteTimes(end) == 1)
                        Transitions = [1; Transitions; length(NoteTimes)];
                    else
                        Transitions = [1; Transitions];
                    end
                else
                    if (NoteTimes(end) == 1)
                        Transitions = [Transitions; length(NoteTimes)];
                    else
                        Transitions = Transitions;
                    end
                end

                ActualTransitions = reshape(Transitions, 2, length(Transitions)/2);
                ShortIntervals = find((ActualTransitions(1,2:end) - ActualTransitions(2,1:(end-1))) < (0.002 * Fs));
                for j = 1:length(ShortIntervals),
                    Transitions(2*(ShortIntervals(j))) = -100;
                    Transitions(2*(ShortIntervals(j)) + 1) = -100;
                end
                Transitions = Transitions(find(Transitions > 0));
                ActualTransitions = reshape(Transitions,2,length(Transitions)/2);
                Onsets = ActualTransitions(1,:)/Fs;
                Offsets = ActualTransitions(2,:)/Fs;

                InterSyllableIntervals = Onsets(2:end) - Offsets(1:(end - 1));
                if (size(InterSyllableIntervals,1) == 1);
                    InterSyllableIntervals = InterSyllableIntervals';
                end


                BoutTransitions = find(InterSyllableIntervals > 0.4);
                BoutTransitions = [0; BoutTransitions; length(Onsets)];
                for BoutNo = 2:length(BoutTransitions),
                    Bouts(BoutIndex).SyllableOnsets = Onsets((BoutTransitions(BoutNo - 1) + 1):BoutTransitions(BoutNo));
                    Bouts(BoutIndex).SyllableOffsets = Offsets((BoutTransitions(BoutNo - 1) + 1):BoutTransitions(BoutNo));
                    Bouts(BoutIndex).Duration = Offsets(BoutTransitions(BoutNo)) - Onsets((BoutTransitions(BoutNo - 1) + 1));
                    Bouts(BoutIndex).FileName = DataFiles(i).name;
                    Bouts(BoutIndex).Onset = Onsets((BoutTransitions(BoutNo - 1) + 1));
                    Bouts(BoutIndex).Offset = Offsets(BoutTransitions(BoutNo));
                    Bouts(BoutIndex).NoofSyllables = length(Bouts(BoutIndex).SyllableOnsets);
                    BoutIndex = BoutIndex + 1;
                end
                disp(['Finished analysing ',DataFiles(i).name]);
            end
        catch
            disp(['Problem analysing ',DataFiles(i).name]);
        end
    end
end

        
        

