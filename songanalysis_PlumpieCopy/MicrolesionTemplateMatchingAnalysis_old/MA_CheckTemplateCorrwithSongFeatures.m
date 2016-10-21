function [] = MA_CheckTemplateCorrwithSongFeatures(ParameterFile, OutputDir, SAPFeatureValueFile, Syllables, TemplateFileName, PreOrPost, Context, Day, TemplateFileSAPFeatureValueFile, Colour)

FileSep = filesep;
Colours = 'rgbcmk';

%====== Load and extract parameters =======================================
disp('Extracting parameters ...');
Parameters = MA_ParseParametersFile(ParameterFile);
%==========================================================================

%======Now extract all the song file names=================================
disp('Extracting song file names ...');
% First for the pre-treatment days
for i = 1:Parameters.NoPreDays,
    disp(['Extracting song files for pre-treatment day #', num2str(i), ' ...']);
    Parameters.PreDirSongFileNames{i} = MA_ExtractSongFileNames(Parameters.PreDirSongFileList{i});
    Parameters.PreUnDirSongFileNames{i} = MA_ExtractSongFileNames(Parameters.PreUnDirSongFileList{i});
end

% Next for the post-treatment days
for i = 1:Parameters.NoPostDays,
    disp(['Extracting song files for post-treatment day #', num2str(i), ' ...']);    
    Parameters.PostDirSongFileNames{i} = MA_ExtractSongFileNames(Parameters.PostDirSongFileList{i});
    Parameters.PostUnDirSongFileNames{i} = MA_ExtractSongFileNames(Parameters.PostUnDirSongFileList{i});
end
%==========================================================================

%====================== Loading templates =================================
disp('Loading motif template ...');
Parameters.SyllableTemplate = load(Parameters.SyllableTemplateFileName);

%==========================================================================

%==========================================================================
[SyllableTemplateDir, SyllableTemplateFileName, SyllableTemplateExt] = fileparts(Parameters.SyllableTemplateFileName);
TemplateMatchOutputDir = fullfile(OutputDir, [SyllableTemplateFileName, SyllableTemplateExt, '.TemplateMatchResults']);
SavedDataFile = fullfile(TemplateMatchOutputDir, [SyllableTemplateFileName, SyllableTemplateExt, '.TemplateMatchResults.SavedData.mat']);

if (~exist(SavedDataFile, 'file'))
    [Parameters] = MA_ExtractTemplateMatchingData(Parameters, OutputDir);
end
load(SavedDataFile);
%==========================================================================
% Now get all the syllables in the template
for SyllTemp = 1:length(Parameters.SyllableTemplate.SyllableTemplates),
    SyllLabel(SyllTemp) = Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label;
end

% Load SAP Feature file
load(SAPFeatureValueFile);
SAPFeatures = handles.ASSL;

load(TemplateFileSAPFeatureValueFile);
TemplateSAPFeatures = handles.ASSL;

Features = [{'Duration'} {'MeanFrequency'} {'Entropy'} {'LogAmplitude'} {'PitchGoodness'} {'FrequencyModulation'} {'EntropyVariance'} {'AmplitudeModulation'} {'FundamentalFrequency'}];
% Now get the template values
TemplateFile = find(cellfun(@length, strfind(TemplateSAPFeatures.FileName, TemplateFileName)));
for i = 1:length(Syllables),
    Indices = find(TemplateSAPFeatures.SyllLabels{TemplateFile} == Syllables(i));
    TemplateSyll_Index(i) = Indices(min([3, length(Indices)]));
    % Get all template Features
    for j = 1:length(Features),
        eval(['Template{i}.', Features{j}, ' = ', 'TemplateSAPFeatures.', Features{j}, '{TemplateFile}(TemplateSyll_Index(i));']);
    end
end

% Now I need to get the match values for each of these syllables
for i = 1:length(Syllables),
    Index = 1;
    Templates_SyllIndex = find(SyllLabel == Syllables(i));
    for j = 1:length(SAPFeatures.SyllOnsets),
        % First check if there is this file in the template match list
        FileIndex = find(cellfun(@length, strfind(eval(['Parameters.', PreOrPost, Context, 'SongFileNames{', num2str(Day), '}']), SAPFeatures.FileName{j})));
        if (~isempty(FileIndex))
            SyllableIndices = find(SAPFeatures.SyllLabels{j} == Syllables(i));
            TemplateMatchValues = eval(['Parameters.Syllable', PreOrPost, Context, 'Results{', num2str(Day), '}{', num2str(Templates_SyllIndex), '}{', num2str(FileIndex), '}']);
            if (~isempty(TemplateMatchValues))
                for k = 1:length(SyllableIndices),
                    if ((j == TemplateFile) && (SyllableIndices(k) == TemplateSyll_Index(i)))
                        continue;
                    end
                    [ClosestMatch_MatchVal, ClosestMatch_Index] = min(abs(TemplateMatchValues(:,2) - SAPFeatures.SyllOnsets{j}(SyllableIndices(k))/1000));
                    if (ClosestMatch_MatchVal < 0.01)
                        MatchVal{i}(Index) = TemplateMatchValues(ClosestMatch_Index,1);
                        for FeatureNo = 1:length(Features),
                            eval(['Syll_FeatValues{i}.', Features{FeatureNo}, '(Index) = ', 'SAPFeatures.', Features{FeatureNo}, '{j}(SyllableIndices(', num2str(k), '));']);
                        end
                        Index = Index + 1;
                    end
                end
            end
        end
    end
end

for i = 1:length(Syllables),
    figure(i);
    set(gcf, 'Position', [427 30 800 700]);
    set(gcf, 'Color', 'w');
    for j = 1:6,
        subplot(3,2,j);
        hold on;
        eval(['plot(','[Syll_FeatValues{i}.', Features{j}, ']/Template{i}.', Features{j}, ', MatchVal{i},', '[''', Colour, 'o''', '])']);
        [r, p] = eval(['corrcoef(','[Syll_FeatValues{i}.', Features{j}, ']/Template{i}.', Features{j}, ', MatchVal{i})']);
        title([Features{j}, ': r = ', num2str(r(1,2)), '; p = ', num2str(p(1,2))], 'FontSize', 12, 'FontWeight', 'bold');
        set(gca, 'FontSize', 12, 'FontWeight', 'bold');
        
        if (mod(j,2) == 1)
            ylabel('Match value', 'FontSize', 12, 'FontWeight', 'bold');
        end
        if (j > 4)
            xlabel('Syllable Feature value / Template Feature value', 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
    annotation('textbox', [0.1 0.97 0.8 0.02], 'String', ['Syllable ', Syllables(i)], 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'LineStyle', 'none')
end
        
disp('Finished');

