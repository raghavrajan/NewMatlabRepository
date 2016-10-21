function [] = MakeSpecificSyllableTemplatesFromFileSAPFeatures(DirectoryName, NoteDir, SongFile, FileType, OutputDir, TemplateIncrement, Sequence, TimeStretch, FreqStretch)

PresentDir = pwd;
cd(DirectoryName);
disp(SongFile);

if (strfind(FileType, 'okrank'))
    [Song, Fs] = ReadOKrankData(DirectoryName, SongFile, 1);
else
    if (strfind(FileType, 'obs'))
        if ispc
            [Song, Fs] = soundin_copy([DirectoryName, '\'], SongFile, 'obs0r');
        else
            [Song, Fs] = soundin_copy([DirectoryName, '/'], SongFile, 'obs0r');
        end
        Song = Song*5/32768;
    else
        if (strfind(FileType, 'wav'))
            cd(DirectoryName);
            [Song, Fs] = wavread(SongFile);
            cd(PresentDir);
        end
    end
end

[m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude ,m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight] = deriv(Song, Fs);

SAPFeatures(1,:) = m_AM(:)';
SAPFeatures(2,:) = m_FM(:)';
SAPFeatures(3,:) = m_Entropy(:)';
SAPFeatures(4,:) = m_amplitude(:)';
SAPFeatures(5,:) = m_Freq(:)';
SAPFeatures(6,:) = m_PitchGoodness(:)';

cd(NoteDir);
load([SongFile, '.not.mat']);
Matches = strfind(labels, Sequence);
for i = 1:length(Matches),
    StartIndex = Matches(i);
    EndIndex = Matches(i) + length(Sequence) - 1;
    
    MotifTemplate = MakeTemplatesMicrolesionSpectralMatchAnalysisSAPFeatures(SAPFeatures, size(SAPFeatures,2)/(length(Song)/Fs), [(onsets(StartIndex)/1000 - 0.01) (offsets(EndIndex)/1000 + 0.01)], Sequence, TimeStretch, FreqStretch);
    
    OutputFileName = [SongFile, '.', Sequence, '.', num2str(i + TemplateIncrement), '.SAPFeatures.template.mat'];
    save([OutputDir, '/', OutputFileName], 'MotifTemplate');
    disp(['Saved templates to ', OutputFileName, ' in directory ', OutputDir]);
end

