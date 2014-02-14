function [SAPFeatures] = CalculateSAPFeatures(SoundFile,Fs)

[File, Fs]=wavread(SoundFile);

%[m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude , m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ]=segment_data(File,Fs);
[m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude , m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ]=deriv(File, Fs);

SAPFeatures.AM = m_AM;
SAPFeatures.Entropy = m_Entropy;
SAPFeatures.Freq = m_Freq;
SAPFeatures.Pitch = Pitch_chose;
SAPFeatures.FM = m_FM;
SAPFeatures.Amplitude = m_amplitude;
SAPFeatures.PitchGoodness = m_PitchGoodness;
SAPFeatures.Time = linspace(0,(length(File)/Fs),length(m_AM));
SAPFeatures.Time = SAPFeatures.Time';

OutputFile = [SoundFile,'.SAPFeatures.mat'];

%save(OutputFile,'SAPFeatures');

% disp('Finished calculating SAP Features');