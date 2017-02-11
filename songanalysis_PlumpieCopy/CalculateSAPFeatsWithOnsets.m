function [Feats] = CalculateSAPFeatsWithOnsets(Song, Time, Fs, Onsets, Offsets)

[m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude ,m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ]=deriv(Song, Fs);
T = linspace(Time(1), Time(end), length(m_Entropy));

for i = 1:length(Onsets),
    SyllNo = i;
    StartIndex = find(T < Onsets(i), 1, 'last');
    EndIndex = find(T > Offsets(i), 1, 'first');
    Feats(SyllNo, 1) = Offsets(i) - Onsets(i); % Duration 
    Feats(SyllNo, 2) = mean(m_amplitude(StartIndex:EndIndex)); % Amplitude
    Feats(SyllNo, 3) = mean(m_Entropy(StartIndex:EndIndex)); % Entropy 
    Feats(SyllNo, 4) = mean(m_Freq(StartIndex:EndIndex));
    Feats(SyllNo, 5) = mean(m_AM(StartIndex:EndIndex));
    Feats(SyllNo, 6) = mean(m_PitchGoodness(StartIndex:EndIndex));
    Feats(SyllNo, 7) = mean(m_FM(StartIndex:EndIndex));
    Feats(SyllNo, 8) = var(m_Entropy(StartIndex:EndIndex));
end        
