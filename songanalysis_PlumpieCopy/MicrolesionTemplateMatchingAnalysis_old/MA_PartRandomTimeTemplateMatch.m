function [Bout] = MA_PartRandomTimeTemplateMatch(DataDir, SongFile, FileType, MotifTemplate, Label, OutputDir, TemplateType, ShufflePercentage)

PresentDir =pwd;

TimeNow = datestr(now, 'mmddyyHHMMSS');

FileSep = filesep;

if (DataDir(end) ~= FileSep)
    DataDir(end+1) = FileSep;
end

FFTWinSize = MotifTemplate.MotifTemplate(1).FFTWinSize;
FFTWinOverlap = MotifTemplate.MotifTemplate(1).FFTWinOverlap;

% Check to see if analyzed file already exists : if it does then go to
% next file

OutputFileName = [SongFile, '.', Label, '.*.TempMatch.mat'];
cd(OutputDir);
if (~isempty(dir(OutputFileName)))
    return;
end

[Song, Fs] = MA_ReadSongFile(DataDir, SongFile, FileType);

if (isempty(Song))
    return;
end
    
try
    SongTime = (1:1:length(Song))/Fs;

    if (strfind(TemplateType, 'Spectrogram'))
        [P, F, S1, T] = CalculateMultiTaperSpectrogram(Song, Fs, FFTWinSize, FFTWinOverlap, 1.5);

        Freq1 = find((F >= 860) & (F <= 8600));
        S = log10(abs(S1(Freq1,:)));
    else
        [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude ,m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight] = deriv(Song, Fs);

        SAPFeatures(1,:) = m_AM(:)';
        SAPFeatures(2,:) = m_FM(:)';
        SAPFeatures(3,:) = m_Entropy(:)';
        SAPFeatures(4,:) = m_amplitude(:)';
        SAPFeatures(5,:) = m_Freq(:)';
        SAPFeatures(6,:) = m_PitchGoodness(:)';

        S = SAPFeatures;
        T = linspace(SongTime(1), SongTime(end), size(S, 2));
    end

    for Repetition = 1:2,
        TempS = zeros(size(S));
        RandomColIndices = randperm(size(S,2));
        TempS = S(:, RandomColIndices);
        for i = 1:length(MotifTemplate.MotifTemplate),
            WMotif1 = MotifTemplate.MotifTemplate(i).MotifTemplate;
            TempMeanSTD = CalculateMeanSTDforSpectralMatch(TempS(1:size(S,1)*size(S,2)), size(WMotif1,1)*size(WMotif1,2), (size(S,2) - size(WMotif1,2) + 1), size(WMotif1,1));

            WinMean = TempMeanSTD(1:length(TempMeanSTD)/2);
            WinSTD = TempMeanSTD((length(TempMeanSTD)/2 + 1):end);

            [Match] = CalTemplateMatch(WMotif1, TempS, WinMean, WinSTD);
            Match = Match*size(WMotif1,1)*size(WMotif1,2);
            Bout.BoutSeqMatch{i} = Match;
        end
        Bout.T = T;

        clear Match;
        for MatchNo = 1:length(Bout.BoutSeqMatch),
            Match(MatchNo,:) = Bout.BoutSeqMatch{MatchNo}(1:length(Bout.BoutSeqMatch{end}));
        end
        Match = max(Match);
        Bout.MaxBoutSeqMatch = Match;
        clear Match;
        [MaxVal, MaxInd] = max(Bout.MaxBoutSeqMatch);
        Bout.MaxBoutSeqMatchVal = [Bout.T(MaxInd) MaxVal];
        Bout.FileName = SongFile;
        Bout.FileLength = SongTime(end);
        Bout.BoutSeqMatch = [];

        cd(OutputDir);
        OutputFileName = [SongFile, '.', Label, '.', num2str(Repetition), '.TempMatch.mat'];
        save(OutputFileName, 'Bout');

        cd(PresentDir);
        clear Bout;
        clear onsets offsets Bouts BoutOnsets BoutOffsets Temp;
    end
catch
    disp(['Could not analyze ', SongFile]);
end

