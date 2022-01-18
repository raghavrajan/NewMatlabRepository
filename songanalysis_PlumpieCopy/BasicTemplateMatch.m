function [Match] = BasicTemplateMatch(RawDataDirectory, SongFile, FileType, TemplateFile)

PresentDir = pwd;

% First load up templates
Templates{1} = load(TemplateFile);


% Now display templates for each of the syllables in a figure;
figure;
for j = 1:length(Templates{1}.SyllableTemplates),
    subplot(2, ceil(length(Templates{1}.SyllableTemplates)/2),j);
    contourf(Templates{1}.SyllableTemplates{j}{1}.MotifTemplate(1).MotifTemplate)
    colorbar;
end

cd(PresentDir);

% Load up the file that has to be matched and calculate the appropriate
% spectrogram with or without raw data normalization
[RawData, Fs] = GetData(RawDataDirectory, SongFile, FileType, 0);

% Now match them one by one to the data file specified
for i = 1:length(Templates),
    fprintf('Template #%d\n', i);
    
    % Get FFTWinSize, FFTWinOverlap, etc.
    FFTWinSize = Templates{i}.SyllableTemplates{1}{1}.MotifTemplate(1).FFTWinSize;
    FFTWinOverlap = Templates{i}.SyllableTemplates{1}{1}.MotifTemplate(1).FFTWinOverlap;
    
    % Now calculate spectrogram
%     if (strfind(Templates{i}.SyllableTemplates{1}{1}.MotifTemplate(1).SpectrogramType, 'multitaper'))
%         [P, F, S1, T] = CalculateMultiTaperSpectrogram(Data, Fs, FFTWinSize, FFTWinOverlap, 1.5);
%     else
    [S1, F, T, P] = spectrogram(RawData, hamming(round(FFTWinSize*Fs/1000)), round(FFTWinOverlap*Fs/1000), round(FFTWinSize*Fs/1000), Fs);
    %end
    
    Freq1 = find((F >= 860) & (F <= 8600));
    Spect = log10(abs(S1(Freq1,:)));
    
    % Now go through and match with each of the different syllable types
    for j = 1:length(Templates{i}.SyllableTemplates),
        fprintf('Syllable #%d >> ', j);
        % Now do this for each of the time and freq stretch variants of the
        % template
        
        Index = 1;
        for k = 1:length(Templates{i}.SyllableTemplates{j}{1}.MotifTemplate),
            Motif = Templates{i}.SyllableTemplates{j}{1}.MotifTemplate(k).MotifTemplate;
            
            TempMeanSTD = CalculateMeanSTDforSpectralMatch(Spect(:), size(Motif,1)*size(Motif,2),  (size(Spect,2) - size(Motif,2) + 1), size(Motif,1));
            WinMean = TempMeanSTD(1:length(TempMeanSTD)/2);
            WinSTD = TempMeanSTD((length(TempMeanSTD)/2) + 1:end);
            ZeroMean = zeros(size(WinMean));
            OneSTD = ones(size(WinSTD));
                        
            TempMatch{Index} = CalTemplateMatch(Motif, Spect, WinMean, WinSTD);
            TempMatch{Index} = TempMatch{Index}*size(Motif,1)*size(Motif,2);
            
            Index = Index + 1;
        end
        
        MinMatchLen = min(cellfun(@length, TempMatch));
        for IndexNum = 1:length(TempMatch),
            Match{i}{j}(IndexNum,:) = TempMatch{IndexNum}(1:MinMatchLen);
        end
        Match{i}{j} = max(Match{i}{j});
    end
    fprintf('\n');
end

disp('Finished');


