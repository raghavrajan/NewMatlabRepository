function [LogAmplitude] = ASSLCalculateLogAmplitudeKaoWithOnsets(Song, Time, Fs, Onsets, Offsets)

if (~isempty(Onsets))
    [SongLogAmplitude] = ASSLCalculateLogAmplitudeKao(Song, Fs, Time, [], []);
    T = Time;
    T = T(:);
    Time = Time(:);
    SongLogAmplitude = SongLogAmplitude(:);
end

if (isempty(Onsets))
    LogAmplitude = []; % Amplitude
end

for i = 1:length(Onsets),
    SyllNo = i;
    StartIndex = find(T <= Onsets(i), 1, 'last');
    if (isempty(StartIndex))
        StartIndex = 1;
    end
    
    EndIndex = find(T >= Offsets(i), 1, 'first');
    if (isempty(EndIndex))
        EndIndex = length(T);
    end
    
%    LogAmplitude(SyllNo) = polyarea(Time(StartIndex:EndIndex), SongLogAmplitude(StartIndex:EndIndex)); % Amplitude
    LogAmplitude(SyllNo) = mean(SongLogAmplitude(StartIndex:EndIndex)); % Amplitude
end        
disp('Finished');