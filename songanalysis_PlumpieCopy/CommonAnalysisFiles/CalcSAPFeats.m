function [SAPFeatsMatrix, SAPFeat_FieldNames] = CalcSAPFeats(BirdParameters)

% This function is to calculate SAP Features for syllables one at a time.
% The idea is to take the sound waveform for each syllable with 50ms
% padding on either side and then use the deriv function from SAM to
% calculate SAP features.

% For continuous data, also, the same function can be used - only when a
% syllable is spread over two files, I have to get the sound waveform from
% both files and then use deriv function

% =========================================================================

tic
SAPFeatsMatrix = [];
Padding = 0.2; % A variable that decides how much extra raw data to take before and after each syllable in seconds
fprintf('\n>');
UniqueFileIndices = unique(BirdParameters.SyllableData(:,2));
for i = 1:length(UniqueFileIndices),
    if (mod(i, 5) == 0)
        fprintf('%d>>', i);
    end
    FileIndices = find(BirdParameters.SyllableData(:,2) == UniqueFileIndices(i));
    
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
    
    Time = (1:1:length(Song))/Fs;
    
    % Now for all syllables that don't have enough padding ahead of the
    % syllable or at the end of the syllable, I basically make all SAP
    % Feats for that syllable as NaNs and then throw that syllable out and
    % then calculate sap features for the rest of the syllables
    
    OnsetsTobeRemoved = find(SyllOnsetTimes <= (Padding*1000)); % have to multiple Padding by 1000 as OnsetTimes are in ms and Padding is in seconds
    OffsetsTobeRemoved = find(SyllOffsetTimes >= ((length(Song)*1000/Fs) - Padding*1000));
    
    SyllOnsetTimes([OnsetsTobeRemoved(:); OffsetsTobeRemoved]) = [];
    SyllOffsetTimes([OnsetsTobeRemoved(:); OffsetsTobeRemoved]) = [];
    
% This part of the code calculates the sap features for the entire file
% with all syllable onsets - took 176.135522 seconds for o19b78-17032011

%     SongForSAPCalculations = Song((round(SyllOnsetTimes(1)*Fs/1000) - Padding*1000):(round(SyllOffsetTimes(end)*Fs/1000) + Padding*1000));
%     TimeForSAPCalculations = Time((round(SyllOnsetTimes(1)*Fs/1000) - Padding*1000):(round(SyllOffsetTimes(end)*Fs/1000) + Padding*1000));
%     
%     [Feats, RawFeats, FeatsFs] = ASSLCalculateSAPFeatsWithOnsets(SongForSAPCalculations, TimeForSAPCalculations, Fs, SyllOnsetTimes/1000, SyllOffsetTimes/1000);
%     
%     if (~exist('SAPFeat_FieldNames', 'var'))
%         SAPFeat_FieldNames = fieldnames(Feats);
%     end
%     
%     TempSAPFeatsMatrix = [];
%     for k = 1:length(SAPFeat_FieldNames),
%         SAPFeats.(SAPFeat_FieldNames{k}) = Feats.(SAPFeat_FieldNames{k})(:);
%         TempSAPFeatsMatrix = [TempSAPFeatsMatrix Feats.(SAPFeat_FieldNames{k})(:)];
%     end
% 
%     SAPFeatsMatrix = [SAPFeatsMatrix; [ones(length(OnsetsTobeRemoved), length(SAPFeat_FieldNames))*NaN; TempSAPFeatsMatrix; ones(length(OffsetsTobeRemoved), length(SAPFeat_FieldNames))*NaN]]; 

% This part of the code calculates the sap features one syllable at a time
% much faster - this took only 37.8917 seconds for o19b78-17032011

    if (~exist('SAPFeat_FieldNames', 'var'))
        SAPFeatsMatrix = [SAPFeatsMatrix; [ones(length(OnsetsTobeRemoved), 9)*NaN]];
    else
        SAPFeatsMatrix = [SAPFeatsMatrix; [ones(length(OnsetsTobeRemoved), length(SAPFeat_FieldNames))*NaN]];
    end
    for j = 1:length(SyllOnsetTimes),
        SongForSAPCalculations = Song((round(SyllOnsetTimes(j)*Fs/1000) - Padding*1000):(round(SyllOffsetTimes(j)*Fs/1000) + Padding*1000));
        TimeForSAPCalculations = Time((round(SyllOnsetTimes(j)*Fs/1000) - Padding*1000):(round(SyllOffsetTimes(j)*Fs/1000) + Padding*1000));
    
        % Again a dirty hack to get this to work for red02yellow06
        try
            [Feats, RawFeats, FeatsFs] = ASSLCalculateSAPFeatsWithOnsets(SongForSAPCalculations, TimeForSAPCalculations, Fs, SyllOnsetTimes(j)/1000, SyllOffsetTimes(j)/1000);
        catch
            Feats = NaN;
            RawFeats = NaN;
            FeatsFs = NaN;
        end
        
        if (isstruct(Feats))
            if (~exist('SAPFeat_FieldNames', 'var'))
                SAPFeat_FieldNames = fieldnames(Feats);
            end

            TempSAPFeatsMatrix = [];
            for k = 1:length(SAPFeat_FieldNames),
                SAPFeats.(SAPFeat_FieldNames{k}) = Feats.(SAPFeat_FieldNames{k})(:);
                TempSAPFeatsMatrix = [TempSAPFeatsMatrix Feats.(SAPFeat_FieldNames{k})(:)];
            end
            SAPFeatsMatrix = [SAPFeatsMatrix; [TempSAPFeatsMatrix]];
        else
            SAPFeatsMatrix = [SAPFeatsMatrix; ones(1,9)*NaN];
        end            
    end
    if (~exist('SAPFeat_FieldNames', 'var'))
        SAPFeatsMatrix = [SAPFeatsMatrix; [ones(length(OffsetsTobeRemoved), 9)*NaN]]; 
    else
        SAPFeatsMatrix = [SAPFeatsMatrix; [ones(length(OffsetsTobeRemoved), length(SAPFeat_FieldNames))*NaN]]; 
    end
end
fprintf('\n');