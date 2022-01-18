function [Correlation] = CalculateMotifTemporalSimilarity(DirectoryName, FileList, NoteFileDir, FileType, Motif)
 
% =========================================================================
% Script for calculating temporal similarity of the amplitude envelope of a
% motif. 
% Calculates the correlation as follows:
% Correlation =    (Motif1 - mean(Motif1)) * (Motif2 - mean(Motif2))'
%               ---------------------------------------------------------
%               norm(Motif1 - mean(Motif1)) * norm(Motif2 - mean(Motif2))
% Inputs:
%       1. DirectoryName - directory where raw data files are stored
%       2. FileList - list of files. Here the first motif in the first file
%       will be used as the template and will be compared to all motifs in
%       the remaining files
%       3. NoteFileDir - directory where notefiles are kept. Doesn't have
%       to be the full path, but this will be considered relative to
%       DirectoryName
%       4. FileType - type of raw data file
%       5. Motif - this gives the labels for the motif that has to be
%       correlated

% Outputs:
%       1. Correlation - this gives the values of the correlations for the
%       template with each of the other motifs
%
% In each case, the template is shrunk by upto 5% and extended by upto 10%
% and the correlation is calculated for each of these versions of the
% template and the highest correlation value is considered for each motif.
% THis is to account for the fact that different motifs are of slightly
% different lengths.
%
% =========================================================================

Correlation = [];
PresentDir = pwd;
cd(DirectoryName);

% First get all the files
Fid = fopen(FileList, 'r');
Files = textscan(Fid, '%s', 'DeLimiter', '\n');
Files = Files{1};
fclose(Fid);

% Now load up the note files
for i = 1:length(Files),
    NoteInfo{i} = load(fullfile(DirectoryName, NoteFileDir, [Files{i}, '.not.mat']));
end

% Now for each of the files, calculate log amplitude (aronov fee) and then
% store this amplitude envelope in a vector

FHigh = 300;
FLow = 8000;
FFTWinSize = 5; % in ms
StretchValues = -5:1:10; % in %

fprintf('Total of %d files\n', length(Files));
fprintf('\n');
Index = 1;
MatchIndex = 1;
for i = 1:length(Files),
    fprintf('%d > ', i);
    [RawData, Fs] = GetData(DirectoryName, Files{i}, FileType, 0);
    FiltData = bandpass(RawData, Fs, FHigh, FLow);
    
    LogAmplitudeAronovFee = ASSLCalculateLogAmplitudeAronovFee(RawData, Fs, [], FFTWinSize, [], FHigh, FLow);
    
    Matches = strfind(NoteInfo{i}.labels, Motif);
    
    if (i == 1)
        if (isempty(Matches))
            disp(['No motifs (', Motif, ') found in the template file']);
            return;
        end
        
        j = Matches(1);
        TemplateAmpWF{1} = LogAmplitudeAronovFee(round(NoteInfo{i}.onsets(j) * Fs/1000): round(NoteInfo{i}.offsets(j+3) * Fs/1000));
        TemplateAmpTime = linspace(1/Fs, length(TemplateAmpWF{1})/Fs, length(TemplateAmpWF{1}));
        
        % Now make all the stretched versions of the template
        Index = 1;
        for Stretch = StretchValues(:)',
            StretchedTime = linspace(TemplateAmpTime(1), TemplateAmpTime(end), length(TemplateAmpTime) * (1 + Stretch/100));
            StretchedTemplate{Index} = interp1(TemplateAmpTime, TemplateAmpWF{1}, StretchedTime);
            Index = Index + 1;
        end
    else
        for j = Matches(:)',
            MotifAmpWF = LogAmplitudeAronovFee(round(NoteInfo{i}.onsets(j) * Fs/1000): round(NoteInfo{i}.offsets(j+3) * Fs/1000));
            TempCorr = [];
            for Stretch = 1:length(StretchedTemplate),
                Template = StretchedTemplate{Stretch};
                if (length(Template) <= length(MotifAmpWF))
                    MotifAmpWF = MotifAmpWF(1:length(Template));
                else
                    Template = Template(1:length(MotifAmpWF));
                end
                TempCorr(end+1) = ((Template - mean(Template)) * (MotifAmpWF - mean(MotifAmpWF))')/(norm(Template - mean(Template)) * norm(MotifAmpWF - mean(MotifAmpWF)));
            end
            Correlation(MatchIndex) = max(TempCorr);
            MatchIndex = MatchIndex + 1;
        end
    end
end

cd(PresentDir);
disp('Finished calculating temporal amplitude similarity');