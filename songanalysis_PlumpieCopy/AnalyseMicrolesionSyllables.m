function [Feats, UnidentifiedSylls, UniqueLabels] = AnalyseMicrolesionSyllables(Dir, Date, Condition, RawDataDir, FileType)

PresentDir = pwd;

cd(Dir);
MaxMatchFiles = dir('*maxMatch*.txt');

for i = 1:length(MaxMatchFiles),
    if (~isempty(strfind(MaxMatchFiles(i).name, [Date, '_', Condition])))
        Temp = textread(MaxMatchFiles(i).name, '%s', 'delimiter', '\n');
    end
end

MaxFiles = min(25, length(Temp));

% Load note files
for i = 1:MaxFiles,
    Notes{i} = load([Temp{i}, '.not.mat']);
    Notes{i}.onsets = Notes{i}.onsets/1000;
    Notes{i}.offsets = Notes{i}.offsets/1000;
end

% Calculate SAP Features for all these files
cd(RawDataDir);
for i = 1:MaxFiles,
    SongFile = Temp{i};
    disp(SongFile);
    if (strfind(FileType, 'okrank'))
        [Song, Fs] = ReadOKrankData(RawDataDir, SongFile, 1);
    else
        if (strfind(FileType, 'wav'))
            [Song, Fs] = wavread(SongFile);
        else
            if (strfind(FileType, 'obs'))
                channel_string = strcat('obs',num2str(0),'r');
                [Song, Fs] = soundin_copy([RawDataDir, '/'], SongFile, channel_string);

                % Convert to V - 5V on the data acquisition is 32768
                Song = Song * 5/32768;
            end
        end
    end
    Time = (0:1:length(Song)-1)/Fs;
    [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude , m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ]=deriv(Song, Fs);
    FeatureTime = linspace(Time(1), Time(end), length(m_Entropy));
    
    for j = 1:length(Notes{i}.onsets),
        StartIndex = find(FeatureTime <= Notes{i}.onsets(j), 1, 'last');
        EndIndex = find(FeatureTime <= Notes{i}.offsets(j), 1, 'last');
        
        Features.AM{i}(j) = mean(m_AM(StartIndex:EndIndex));
        Features.FM{i}(j) = mean(m_FM(StartIndex:EndIndex));
        Features.Amplitude{i}(j) = mean(m_amplitude(StartIndex:EndIndex));
        Features.Pitch{i}(j) = mean(Pitch_chose(StartIndex:EndIndex));
        Features.PG{i}(j) = mean(m_PitchGoodness(StartIndex:EndIndex));
        Features.Freq{i}(j) = mean(m_Freq(StartIndex:EndIndex));
        Features.Entropy{i}(j) = mean(m_Entropy(StartIndex:EndIndex));
    end
end

Labels = [];

for i = 1:length(Notes),
    Labels = [Labels Notes{i}.labels];
end

UniqueLabels = unique(Labels);

for i = 1:length(UniqueLabels),
    Feats(i).Label = UniqueLabels(i);
    Feats(i).AM = [];
    Feats(i).FM = [];
    Feats(i).Amplitude = [];
    Feats(i).Pitch = [];
    Feats(i).PG = [];
    Feats(i).Freq = [];
    Feats(i).Entropy = [];
    Feats(i).Dur = [];
    
    for j = 1:length(Notes),
        Matches = find(Notes{j}.labels == UniqueLabels(i));
        Feats(i).AM = [Feats(i).AM Features.AM{j}(Matches)];
        Feats(i).FM = [Feats(i).FM Features.FM{j}(Matches)];
        Feats(i).Amplitude = [Feats(i).Amplitude Features.Amplitude{j}(Matches)];
        Feats(i).Pitch = [Feats(i).Pitch Features.Pitch{j}(Matches)];
        Feats(i).PG = [Feats(i).PG Features.PG{j}(Matches)];
        Feats(i).Freq = [Feats(i).Freq Features.Freq{j}(Matches)];
        Feats(i).Entropy = [Feats(i).Entropy Features.Entropy{j}(Matches)];
        Feats(i).Dur = [Feats(i).Dur; (Notes{j}.offsets(Matches) - Notes{j}.onsets(Matches))];
    end
end

UnidentifiedSylls = length(find(Labels == '?'))/length(Labels);

cd(PresentDir);
disp('Finished Analysis');
