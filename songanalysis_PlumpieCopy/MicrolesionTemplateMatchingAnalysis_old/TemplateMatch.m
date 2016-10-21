function [Bout] = TemplateMatch(DirectoryName, FileList, FileType, MotifTemplate, Label, OutputDir, StatusText, StatusAxis, StatusString, Colour, TemplateType)

PresentDir =pwd;

TimeNow = datestr(now, 'mmddyyHHMMSS');

FileSep = filesep;

if (DirectoryName(end) ~= FileSep)
    DirectoryName(end+1) = FileSep;
end

Fid = fopen(FileList, 'r');
Temp = textscan(Fid, '%s', 'delimiter', '\n');
fclose(Fid);

SongFiles = Temp{1};

cd(DirectoryName);

FFTWinSize = MotifTemplate.MotifTemplate(1).FFTWinSize;
FFTWinOverlap = MotifTemplate.MotifTemplate(1).FFTWinOverlap;

%cla(StatusAxis);
%axes(StatusAxis);
%rectangle('Position', [1 0 length(SongFiles)+1 1], 'EdgeColor', 'w', 'FaceColor', 'w');
%axis([1 length(SongFiles)+1 0 1]);
%set(gca, 'Visible', 'off');
%hold on;

for SongFileNo = 1:length(SongFiles),
    SongFile = SongFiles{SongFileNo};

    set(StatusText, 'String', [StatusString, Label, ': ', num2str(SongFileNo), ' of ', num2str(length(SongFiles))]);
    rectangle('Position', [SongFileNo 0 1 1], 'EdgeColor', Colour, 'FaceColor', Colour);
    
    if (isempty(SongFile))
        continue;
    end
    
    Slash = find((SongFile == '/') | (SongFile == '\'));
    if (~isempty(Slash))
        SongFile = SongFile(Slash(end)+1:end);
    end
    disp(SongFile);

    % Check to see if analyzed file already exists : if it does then go to
    % next file
    OutputFileName = [SongFile, '.', Label, '.TempMatch.mat'];
    cd(OutputDir);
    if (exist(OutputFileName, 'file'))
        continue;
    end
    
    cd(DirectoryName);
    
    try
        if (strfind(FileType, 'okrank'))
            [Song, Fs] = ReadOKrankData(DirectoryName, SongFile, 1);
        else
            if (strfind(FileType, 'wav'))
                [Song, Fs] = wavread(SongFile);
            else
                 if (strfind(FileType, 'obs'))
                    channel_string = strcat('obs',num2str(0),'r');
                    [Song, Fs] = soundin_copy(DirectoryName, SongFile, channel_string);
    
                    % Convert to uV - 5V on the data acquisition is 32768
                    Song = Song * 5/32768;
                end
            end
        end
    catch
        continue;
    end
    
    if (isempty(Song))
        continue;
    end

    try
        SongTime = (1:1:length(Song))/Fs;
       
        if (strfind(TemplateType, 'Spectrogram'))
            WinSize = round(FFTWinSize * Fs);
            WinOverlap = round(FFTWinOverlap * WinSize);

            %[S1, F, T, P] = spectrogram(Song, hamming(WinSize), WinOverlap, WinSize, Fs);
            [P, F, S1, T] = CalculateMultiTaperSpectrogram(Song, Fs, 8, 8/2, 1.5);
            
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
        
        for i = 1:length(MotifTemplate.MotifTemplate),
            WMotif = MotifTemplate.MotifTemplate(i).MotifTemplate;
            WinMean = zeros((size(S,2) - size(WMotif,2) + 1), 1);
            WinSTD = zeros((size(S,2) - size(WMotif,2) + 1), 1);

            TempMeanSTD = CalculateMeanSTDforSpectralMatch(S(1:size(S,1)*size(S,2)), size(WMotif,1)*size(WMotif,2), (size(S,2) - size(WMotif,2) + 1), size(WMotif,1));

            WinMean = TempMeanSTD(1:length(TempMeanSTD)/2);
            WinSTD = TempMeanSTD((length(TempMeanSTD)/2 + 1):end);
            [Match] = CalTemplateMatch(WMotif, S, WinMean, WinSTD);
            Match = Match*size(WMotif,1)*size(WMotif,2);
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
   catch
       disp(['Could not analyse ', SongFile]);
    end
   
   cd(OutputDir);
   save(OutputFileName, 'Bout');
   
   cd(DirectoryName);
   clear Bout;
    clear onsets offsets Bouts BoutOnsets BoutOffsets Temp;
end

