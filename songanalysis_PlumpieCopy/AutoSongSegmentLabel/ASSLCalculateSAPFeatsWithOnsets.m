function [Feats, RawFeats] = ASSLCalculateSAPFeatsWithOnsets(Song, Time, Fs, Onsets, Offsets)


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

    RawFeats.LogAmplitude = []; % Amplitude
    RawFeats.Entropy = []; % Entropy 
    RawFeats.MeanFrequency = [];
    RawFeats.AmplitudeModulation = [];
    RawFeats.PitchGoodness = [];
    RawFeats.FrequencyModulation = [];
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
    
    Feats.Duration(SyllNo) = Offsets(i) - Onsets(i); % Duration 
    Feats.LogAmplitude(SyllNo) = mean(m_amplitude(StartIndex:EndIndex)); % Amplitude
    Feats.Entropy(SyllNo) = mean(m_Entropy(StartIndex:EndIndex)); % Entropy 
    Feats.MeanFrequency(SyllNo) = mean(m_Freq(StartIndex:EndIndex));
    Feats.AmplitudeModulation(SyllNo) = mean(m_AM(StartIndex:EndIndex));
    Feats.PitchGoodness(SyllNo) = mean(m_PitchGoodness(StartIndex:EndIndex));
    Feats.FrequencyModulation(SyllNo) = mean(m_FM(StartIndex:EndIndex));
    Feats.EntropyVariance(SyllNo) = var(m_Entropy(StartIndex:EndIndex));
    
    RawFeats.LogAmplitude{SyllNo} = (m_amplitude(StartIndex:EndIndex)); % Amplitude
    RawFeats.Entropy{SyllNo} = (m_Entropy(StartIndex:EndIndex)); % Entropy 
    RawFeats.MeanFrequency{SyllNo} = (m_Freq(StartIndex:EndIndex));
    RawFeats.AmplitudeModulation{SyllNo} = (m_AM(StartIndex:EndIndex));
    RawFeats.PitchGoodness{SyllNo} = (m_PitchGoodness(StartIndex:EndIndex));
    RawFeats.FrequencyModulation{SyllNo} = (m_FM(StartIndex:EndIndex));
end        
