function [Feats, RawFeats, FeatsFs] = ASSLCalculateSAPFeatsWithOnsets_SplitSylls(Song, Time, Fs, Onsets, Offsets)

FeatsToCalculate = [{'Duration'} {'LogAmplitude'} {'Entropy'} {'MeanFrequency'} {'AmplitudeModulation'} {'PitchGoodness'} {'FrequencyModulation'} {'FundamentalFrequency'} {'EntropyVariance'}];
RawFeatsToCalculate = [{'LogAmplitude'} {'Entropy'} {'MeanFrequency'} {'AmplitudeModulation'} {'PitchGoodness'} {'FrequencyModulation'} {'FundamentalFrequency'}];

% if (~isempty(Onsets))
%     
%     % [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude ,m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ]=deriv(Song, Fs);
%     T = linspace(Time(1), Time(end), length(m_Entropy));
%     FeatsFs.SAPFeats_Fs = 1/(T(2) - T(1));
% 
%     % Parameters P for calculating FF for yin
%     % sr for sample rate
%     % maxf0 for max ff
%     % minf0 for min ff
%     P.sr = Fs;
%     P.minf0 = 350; % min ff set to 350 Hz
%     
%     % A dirty fix for some error I was getting - am not sure why the error
%     % was coming so used a try catch routine to bypass this error for the n
%     try
%         FF = yin(Song, P);
%         FF_T = linspace(1/Fs, length(Song)/Fs, length(FF.f0));
%         FF = 2.^FF.f0 * 440;
%         FeatsFs.FF_Fs = 1/(FF_T(2) - FF_T(1));
%     catch
%         FF = NaN;
%         FeatsFs.FF_Fs = NaN;
%     end
% end

if (isempty(Onsets))
    for i = 1:length(FeatsToCalculate),
    	eval(['Feats.', eval(['FeatsToCalculate{', num2str(i), '}']), ' = [];']);
    end
    % Added rms amplitude also in the list of things calculated
    Feats.RMSAmplitude = [];
    % Added a few more features also in the list of things to calculate
    Feats.FirstHalfFM = [];
    Feats.SecondHalfFM = [];
    Feats.FirstHalfPG = [];
    Feats.SecondHalfPG = [];
    Feats.SDMeanFreq = [];
    
    % Feats.Duration = []; % Duration 
    % Feats.LogAmplitude = []; % Amplitude
    % Feats.Entropy = []; % Entropy 
    % Feats.MeanFrequency = [];
    % Feats.AmplitudeModulation = [];
    % Feats.PitchGoodness = [];
    % Feats.FrequencyModulation = [];
    % Feats.EntropyVariance = [];
    % Feats.FundamentalFrequency = [];
    FeatsFs.SAPFeats_Fs = [];
    FeatsFs.FF_Fs = [];
    
    for i = 1:length(RawFeatsToCalculate),
    	eval(['RawFeats.', eval(['RawFeatsToCalculate{', num2str(i), '}']), ' = [];']);
    end
    
%     RawFeats.LogAmplitude = []; % Amplitude
%     RawFeats.Entropy = []; % Entropy 
%     RawFeats.MeanFrequency = [];
%     RawFeats.AmplitudeModulation = [];
%     RawFeats.PitchGoodness = [];
%     RawFeats.FrequencyModulation = [];
%     RawFeats.FundamentalFrequency = [];
    
end

Padding = 0.1; % a variable that decides how much extra raw data to take before and after each syllable in seconds

for i = 1:length(Onsets),
    if (Onsets(i) < Padding)
    	StartIndex = 1;
    else
        % StartIndex = round((Onsets(i) - Padding) * Fs);
        StartIndex = find(Time <= (Onsets(i) - Padding), 1, 'last');
    end
    if (Offsets(i) > (length(Song)/Fs - Padding))
        EndIndex = length(Song);
    else
        % EndIndex = round((Offsets(i) + Padding) * Fs);
        EndIndex = find(Time <= (Offsets(i) + Padding), 1, 'last');
    end

    [TempFeats, TempRawFeats, TempFeatsFs] = ASSLCalculateSAPFeatsWithOnsets(Song(StartIndex:EndIndex), Time(StartIndex:EndIndex), Fs, Onsets(i), Offsets(i));

    FeatsFs = TempFeatsFs;
    
    SyllNo = i;

    for j = 1:length(FeatsToCalculate),
        eval(['Feats.', eval(['FeatsToCalculate{', num2str(j), '}']), '(', num2str(i), ') = TempFeats.', eval(['FeatsToCalculate{', num2str(j), '}']), ';']);
    end
    
    Feats.RMSAmplitude(i) = sqrt(mean(Song(StartIndex:EndIndex).^2));
    
    LenSyll = length(TempRawFeats.FrequencyModulation{1});
    Feats.FirstHalfFM(i) = mean(TempRawFeats.FrequencyModulation{1}(1:round(LenSyll/2)));
    Feats.SecondHalfFM(i) = mean(TempRawFeats.FrequencyModulation{1}(1+round(LenSyll/2):end));
    Feats.FirstHalfPG(i) = mean(TempRawFeats.PitchGoodness{1}(1:round(LenSyll/2)));
    Feats.SecondHalfPG(i) = mean(TempRawFeats.PitchGoodness{1}(1+round(LenSyll/2):end));
    Feats.SDMeanFreq(i) = std(TempRawFeats.MeanFrequency{1});
    
    % Feats.Duration(SyllNo) = Offsets(i) - Onsets(i); % Duration 
    % Feats.LogAmplitude(SyllNo) = mean(m_amplitude(StartIndex:EndIndex)); % Amplitude
    % Feats.Entropy(SyllNo) = mean(m_Entropy(StartIndex:EndIndex)); % Entropy 
    % Feats.MeanFrequency(SyllNo) = mean(m_Freq(StartIndex:EndIndex));
    % Feats.AmplitudeModulation(SyllNo) = mean(m_AM(StartIndex:EndIndex));
    % Feats.PitchGoodness(SyllNo) = mean(m_PitchGoodness(StartIndex:EndIndex));
    % Feats.FrequencyModulation(SyllNo) = mean(m_FM(StartIndex:EndIndex));
    % Feats.EntropyVariance(SyllNo) = var(m_Entropy(StartIndex:EndIndex));
    
    for j = 1:length(RawFeatsToCalculate),
        eval(['RawFeats.', eval(['RawFeatsToCalculate{', num2str(j), '}']), '(', num2str(i), ') = TempRawFeats.', eval(['RawFeatsToCalculate{', num2str(j), '}']), ';']);
    end
    
%     RawFeats.LogAmplitude{SyllNo} = (m_amplitude(StartIndex:EndIndex)); % Amplitude
%     RawFeats.Entropy{SyllNo} = (m_Entropy(StartIndex:EndIndex)); % Entropy 
%     RawFeats.MeanFrequency{SyllNo} = (m_Freq(StartIndex:EndIndex));
%     RawFeats.AmplitudeModulation{SyllNo} = (m_AM(StartIndex:EndIndex));
%     RawFeats.PitchGoodness{SyllNo} = (m_PitchGoodness(StartIndex:EndIndex));
%     RawFeats.FrequencyModulation{SyllNo} = (m_FM(StartIndex:EndIndex));
    
    % if (~isnan(FF))
    %    
    %    StartIndex = find(FF_T <= Onsets(i), 1, 'last');
    %    if (isempty(StartIndex))
    %        StartIndex = 1;
    %    end
    %
    %    EndIndex = find(FF_T >= Offsets(i), 1, 'first');
    %    if (isempty(EndIndex))
    %        EndIndex = length(FF_T);
    %    end
    %
    %    Feats.FundamentalFrequency(SyllNo) = mean(FF(StartIndex:EndIndex));
    %
    %    RawFeats.FundamentalFrequency{SyllNo} = (FF(StartIndex:EndIndex)); % fundamental frequency
    %  else
    %    Feats.FundamentalFrequency(SyllNo) = NaN;
    %    RawFeats.FundamentalFrequency{SyllNo} = NaN;
    % end
end        
