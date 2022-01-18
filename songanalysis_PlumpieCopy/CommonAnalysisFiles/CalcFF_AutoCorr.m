function [FF] = CalcFF_AutoCorr(BirdParameters)

% This function is to calculate FF for syllables one at a time.
% The idea is to take the sound waveform for each syllable with 50ms
% padding on either side and then use the deriv function from SAM to
% calculate SAP features.

% For continuous data, also, the same function can be used - only when a
% syllable is spread over two files, I have to get the sound waveform from
% both files and then use deriv function

% =========================================================================

PlotOption = 0;
FF = [];
Padding = 0.01; % A variable that decides how much extra raw data to take before and after each syllable in seconds
UniqueFileIndices = unique(BirdParameters.SyllableData(:,2));
fprintf('\n>%d>>', length(UniqueFileIndices));
for i = 1:length(UniqueFileIndices), 
    if (mod(i, 5) == 0)
        fprintf('%d>>', i);
    end
    FileIndices = find(BirdParameters.SyllableData(:,2) == UniqueFileIndices(i));
    SyllLabels = char(BirdParameters.SyllableData(FileIndices, 1));
    
    % For each of the syllable labels, find the matching parameter values
    % in the BirdParameters.FFSyllable fields and then extract the start
    % percentage, end percentage, max FF and min FF into a vector the same
    % size as labels
    BirdParameters_FieldNames = fieldnames(BirdParameters);
    FF_Fields = find(cellfun(@length, strfind(BirdParameters_FieldNames, 'FFSyllable')));
    FF_SyllLabels = [];
    
    FF_Index = 1;
    for j = 1:length(FF_Fields),
        FFParameters_String = eval(['BirdParameters.', BirdParameters_FieldNames{FF_Fields(j)}]);
        ColonIndex = find(FFParameters_String == ':');
        if (~isempty(ColonIndex))
            FF_SyllLabels(FF_Index) = FFParameters_String(ColonIndex - 1);
            
            CommaIndices = find(FFParameters_String == ',');
            FF_StartPercentage(FF_Index) = str2double(FFParameters_String(ColonIndex + 1:CommaIndices(1) - 1));
            FF_EndPercentage(FF_Index) = str2double(FFParameters_String(CommaIndices(1) + 1:CommaIndices(2) - 1));
            FF_MaxFF(FF_Index) = str2double(FFParameters_String(CommaIndices(2) + 1:CommaIndices(3) - 1));
            FF_MinFF(FF_Index) = str2double(FFParameters_String(CommaIndices(3) + 1:end));
            FF_Index = FF_Index + 1;
        end
    end
    for j = 1:length(SyllLabels),
        FF_SyllIndex = find(FF_SyllLabels == SyllLabels(j));
        if (isempty(FF_SyllIndex))
            StartPercentage(j) = 0;
            EndPercentage(j) = 100;
            MinFF(j) = 300;
            MaxFF(j) = 2000;
        else
            StartPercentage(j) = FF_StartPercentage(FF_SyllIndex);
            EndPercentage(j) = FF_EndPercentage(FF_SyllIndex);
            MinFF(j) = FF_MinFF(FF_SyllIndex);
            MaxFF(j) = FF_MaxFF(FF_SyllIndex);
        end
    end
    % Check if the offsets for all these syllables is in the same file or not
    if ((mean(BirdParameters.SyllableData(FileIndices, 3)) == UniqueFileIndices(i)) || (UniqueFileIndices(i) == length(BirdParameters.SongFileNames)))
        [Song, Fs] = GetData(BirdParameters.DataDirectory, BirdParameters.SongFileNames{UniqueFileIndices(i)}, BirdParameters.FileType, 0);
        SyllOnsetTimes = BirdParameters.SyllableData(FileIndices,4); % in ms
        SyllOffsetTimes = BirdParameters.SyllableData(FileIndices, 5); % in ms
    else
        [Song1, Fs] = GetData(BirdParameters.DataDirectory, BirdParameters.SongFileNames{UniqueFileIndices(i)}, BirdParameters.FileType, 0);
        [Song2, Fs] = GetData(BirdParameters.DataDirectory, BirdParameters.SongFileNames{UniqueFileIndices(i)+1}, BirdParameters.FileType, 0);
        Song = [Song1(:); Song2(:)];
        SyllOnsetTimes = BirdParameters.SyllableData(FileIndices, 4); % in ms
        % There can never be more than one syllable with an onset in one
        % file and an offset in the next file. So I'm going to take offsets
        % for the first to second last syllable from the first file and
        % offset for the last syllable from the second file
        SyllOffsetTimes = [BirdParameters.SyllableData(FileIndices(1:end-1),5); (BirdParameters.SyllableData(FileIndices(end),5) + (length(Song1)*1000/Fs))]; % in ms
    end
    
    Song = bandpass(Song, Fs, 300, 8000, 'butter');
    Time = (1:1:length(Song))/Fs; % in sec
    
    % Now for all syllables that don't have enough padding ahead of the
    % syllable or at the end of the syllable, I basically make all SAP
    % Feats for that syllable as NaNs and then throw that syllable out and
    % then calculate sap features for the rest of the syllables
    
    OnsetsTobeRemoved = find(SyllOnsetTimes <= (Padding*1000)); % have to multiple Padding by 1000 as OnsetTimes are in ms and Padding is in seconds
    OffsetsTobeRemoved = find(SyllOffsetTimes >= ((length(Song)*1000/Fs) - Padding*1000));
    
    SyllOnsetTimes([OnsetsTobeRemoved(:); OffsetsTobeRemoved]) = [];
    SyllOffsetTimes([OnsetsTobeRemoved(:); OffsetsTobeRemoved]) = [];
    SyllLabels([OnsetsTobeRemoved(:); OffsetsTobeRemoved]) = [];
    StartPercentage([OnsetsTobeRemoved(:); OffsetsTobeRemoved]) = [];
    EndPercentage([OnsetsTobeRemoved(:); OffsetsTobeRemoved]) = [];
    MaxFF([OnsetsTobeRemoved(:); OffsetsTobeRemoved]) = [];
    MinFF([OnsetsTobeRemoved(:); OffsetsTobeRemoved]) = [];
    
    FF = [FF; [ones(length(OnsetsTobeRemoved),1)*NaN]];
    for j = 1:length(SyllOnsetTimes),
        Syllable = Song((round(SyllOnsetTimes(j) * Fs/1000)):(round(SyllOffsetTimes(j) * Fs/1000)));
        NoteSegment = Song((round((SyllOnsetTimes(j) + ((StartPercentage(j)/100) * (SyllOffsetTimes(j) - SyllOnsetTimes(j)))) * Fs/1000)):(round((SyllOnsetTimes(j) + ((EndPercentage(j)/100) * (SyllOffsetTimes(j) - SyllOnsetTimes(j)))) * Fs/1000)));
        if ((~isempty(find(char(FF_SyllLabels) == SyllLabels(j)))) && (~isempty(find(BirdParameters.MotifLabels == SyllLabels(j)))))
            try
                FF(end+1,1) = autocorr_ff_pinterp_plotnotes_WithoutPrompts(Syllable, NoteSegment, Fs, MinFF(j), MaxFF(j), StartPercentage(j), EndPercentage(j), PlotOption);
                CurrentFig = gcf;
                if (PlotOption == 1)
                    set(CurrentFig, 'Position', [153 440 450 650]);
                    %set(gca, 'Position', [0.25 0.25 0.6 0.6]);
                    title([BirdParameters.BirdName, ':', BirdParameters.DataLabel, ':', BirdParameters.Condition, ':', SyllLabels(j), ':', num2str(round(10*FF(end,1))/10), ' Hz:', num2str(round(10*(SyllOffsetTimes(j) - SyllOnsetTimes(j)))/10), ' ms']);
                    if (CurrentFig.Number >= 25)
                        set(CurrentFig, 'Position', [675 440 450 650]);
                        uiwait(CurrentFig);
                        close all;
                        %PlotOption = 0;
                    end
                else
                    close(CurrentFig);
                end
            catch
                CurrentFig = gcf;
                close(CurrentFig);
                if ((~isempty(find(char(FF_SyllLabels) == SyllLabels(j)))) && (~isempty(find(BirdParameters.MotifLabels == SyllLabels(j)))))
                    disp(['Skipped syllable ', SyllLabels(j)]);
                end
                FF(end+1,1) = NaN;
            end
        else
            FF(end+1,1) = NaN;
        end
    end
    FF = [FF; [ones(length(OffsetsTobeRemoved), 1)*NaN]]; 
end
fprintf('\n');