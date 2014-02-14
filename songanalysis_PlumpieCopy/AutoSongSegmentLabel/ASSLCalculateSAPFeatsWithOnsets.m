function [Feats] = ASSLCalculateSAPFeatsWithOnsets(Song, Time, Fs, Onsets, Offsets)


[m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude ,m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ]=deriv(Song, Fs);
T = linspace(Time(1), Time(end), length(m_Entropy));

if (isempty(Onsets))
    Feats.Duration = []; % Duration 
    Feats.LogAmplitude = []; % Amplitude
    Feats.Entropy = []; % Entropy 
    Feats.MeanFrequency = [];
    Feats.AmplitudeModulation = [];
    Feats.PitchGoodness = [];
    Feats.FrequencyModulation = [];
    Feats.EntropyVariance = [];
end
for i = 1:length(Onsets),
    SyllNo = i;
    StartIndex = find(T < Onsets(i), 1, 'last');
    EndIndex = find(T > Offsets(i), 1, 'first');
    Feats.Duration(SyllNo) = Offsets(i) - Onsets(i); % Duration 
    Feats.LogAmplitude(SyllNo) = mean(m_amplitude(StartIndex:EndIndex)); % Amplitude
    Feats.Entropy(SyllNo) = mean(m_Entropy(StartIndex:EndIndex)); % Entropy 
    Feats.MeanFrequency(SyllNo) = mean(m_Freq(StartIndex:EndIndex));
    Feats.AmplitudeModulation(SyllNo) = mean(m_AM(StartIndex:EndIndex));
    Feats.PitchGoodness(SyllNo) = mean(m_PitchGoodness(StartIndex:EndIndex));
    Feats.FrequencyModulation(SyllNo) = mean(m_FM(StartIndex:EndIndex));
    Feats.EntropyVariance(SyllNo) = var(m_Entropy(StartIndex:EndIndex));
end        
