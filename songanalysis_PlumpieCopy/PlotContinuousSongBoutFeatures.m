function [] = PlotContinuousSongBoutFeatures(DirName, FileList, FileType)

PresentDir = pwd;

Fid = fopen(FileList, 'r');
Files = textscan(Fid, '%s', 'DeLimiter', '\n');
Files = Files{1};
fclose(Fid);

figure;

for i = 1:length(Files),
    [Data, Fs] = GetData(DirName, Files{i}, FileType, 0);
    % Data = Data(round(36.3*Fs:39*Fs));
    DataTime = (1:1:length(Data))/Fs;
    
    [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude , m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight] = deriv(Data, Fs);
    
    m_amplitude = zscore(m_amplitude(:));
    m_Entropy = zscore(m_Entropy(:));
    m_Freq = zscore(m_Freq(:));
    m_PitchGoodness = zscore(m_PitchGoodness(:));
    m_FM = zscore(m_FM(:));
    
    [coeff, score, latent] = princomp(zscore([m_amplitude(:), m_Entropy(:), m_Freq(:), m_FM(:), m_PitchGoodness(:)]));
    
    SAPTime = linspace(DataTime(1), DataTime(end), length(m_amplitude));
    
    cd('ASSLNoteFiles');
    NoteInfo = load([Files{i}, '.not.mat']);
    
    Colours = 'rgbcmk';
    
    UniqueSylls = unique(NoteInfo.labels);
    
    for k = 1:length(UniqueSylls)-1,
        Indices = find(NoteInfo.labels == UniqueSylls(k));
        Index = 1;
        for j = [Indices(:)]',
            OnsetTime = find(SAPTime <= (NoteInfo.onsets(j)/1000 - 0.002), 1, 'last');
            OffsetTime = find(SAPTime <= (NoteInfo.offsets(j)/1000 + 0.002), 1, 'last');
            % plot3(score(OnsetTime:OffsetTime,1), score(OnsetTime:OffsetTime,2), score(OnsetTime:OffsetTime,3), [Colours(mod(k-1,length(Colours)) + 1), '.-']);
            if (j == Indices(1))
                plot3(m_amplitude(OnsetTime:OffsetTime), m_Entropy(OnsetTime:OffsetTime), m_FM(OnsetTime:OffsetTime), [Colours(mod(k-1,length(Colours)) + 1), '.-']);
            end
            Features{k}{Index} = [m_amplitude(OnsetTime:OffsetTime) m_Entropy(OnsetTime:OffsetTime) m_Freq(OnsetTime:OffsetTime) m_FM(OnsetTime:OffsetTime) m_PitchGoodness(OnsetTime:OffsetTime)];
            Index = Index + 1;
            hold on;
        end
        
    end
    cd(PresentDir);
    
end
