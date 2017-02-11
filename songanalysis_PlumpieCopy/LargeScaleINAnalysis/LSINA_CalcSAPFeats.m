function [SAPFeats] = LSINA_CalcSAPFeats(BirdParameters)

% This function is to calculate SAP Features for syllables one at a time.
% The idea is to take the sound waveform for each syllable with 50ms
% padding on either side and then use the deriv function from SAM to
% calculate SAP features.

% For continuous data, also, the same function can be used - only when a
% syllable is spread over two files, I have to get the sound waveform from
% both files and then use deriv function

% =========================================================================

SAPFeats = [];
Padding = 0.2; % A variable that decides how much extra raw data to take before and after each syllable in seconds
fprintf('\n>');
UniqueFileIndices = unique(BirdParameters.SyllableData(:,2));
if (BirdParameters.Continuousdata == 0)
    for i = 1:length(UniqueFileIndices),
        if (mod(i, 5) == 0)
            fprintf('%d>>', i);
        end
        FileIndices = find(BirdParameters.SyllableData(:,2) == UniqueFileIndices(i));
        % Check if the offsets for all these syllables is the same or not

        if (BirdParameters.SyllableData(i,2) == BirdParameters.SyllableData(i,3))
            [Song, Fs] = GetData(BirdParameters.DataDirectory, BirdParameters.SongFileNames{BirdParameters.SyllableData(i,2)}, BirdParameters.FileType, 0);

            SongBoutOnsetIndex = find(Time >= (Onsets(1)/1000 - Padding), 1, 'first');
            SongBoutOffsetIndex = find(Time >= (Offsets(end)/1000 + Padding), 1, 'first');

            [Feats, RawFeats, FeatsFs] = ASSLCalculateSAPFeatsWithOnsets(Song(SongBoutOnsetIndex:SongBoutOffsetIndex), Time(SongBoutOnsetIndex:SongBoutOffsetIndex), Fs, Onsets/1000, Offsets/1000);
            if (~exist('SAPFeat_FieldNames', 'var'))
                SAPFeat_FieldNames = fieldnames(Feats);
                for k = 1:length(SAPFeat_FieldNames),
                    SAPFeats.(SAPFeat_FieldNames{k}) = Feats.(SAPFeat_FieldNames{k})(:);
                end
            else
                for k = 1:length(SAPFeat_FieldNames),
                    SAPFeats.(SAPFeat_FieldNames{k}) = [SAPFeats.(SAPFeat_FieldNames{k}); Feats.(SAPFeat_FieldNames{k})(:)];
                end
            end
        end
    end
else
    ContinuousFileTime = [0 cumsum(BirdParameters.FileLen)];
    for j = 1:size(BirdParameters.BoutIndices, 1),
        if (mod(j, 5) == 0)
            fprintf('%d>>', j);
        end
        % Now, first thing is to find out if the onset of the bout and
        % offset of the bout are in the same file or not
        BoutOnsetFile = find(ContinuousFileTime <= (BirdParameters.ActualAllOnsets(BirdParameters.BoutIndices(j,1)) - Padding*1000), 1, 'last');
        BoutOffsetFile = find(ContinuousFileTime <= (BirdParameters.ActualAllOffsets(BirdParameters.BoutIndices(j,2)) + Padding*1000), 1, 'last');
      
        if (BoutOnsetFile == BoutOffsetFile)
            [Song, Fs] = GetData(BirdParameters.DataDirectory, BirdParameters.SongFileNames{BoutOnsetFile}, BirdParameters.FileType, 0);
            Time = (1:1:length(Song))/Fs;
            
            Onsets = BirdParameters.ActualAllOnsets(BirdParameters.BoutIndices(j,1):BirdParameters.BoutIndices(j,2)) - ContinuousFileTime(BoutOnsetFile);
            Offsets = BirdParameters.ActualAllOffsets(BirdParameters.BoutIndices(j,1):BirdParameters.BoutIndices(j,2)) - ContinuousFileTime(BoutOnsetFile);
                
            SongBoutOnsetIndex = find(Time >= (Onsets(1)/1000 - Padding), 1, 'first');
            SongBoutOffsetIndex = find(Time >= (Offsets(end)/1000 + Padding), 1, 'first');
                
            [Feats, RawFeats, FeatsFs] = ASSLCalculateSAPFeatsWithOnsets(Song(SongBoutOnsetIndex:SongBoutOffsetIndex), Time(SongBoutOnsetIndex:SongBoutOffsetIndex), Fs, Onsets/1000, Offsets/1000);
            if (~exist('SAPFeat_FieldNames', 'var'))
                SAPFeat_FieldNames = fieldnames(Feats);
                for k = 1:length(SAPFeat_FieldNames),
                    SAPFeats.(SAPFeat_FieldNames{k}) = Feats.(SAPFeat_FieldNames{k})(:);
                end
            else
                for k = 1:length(SAPFeat_FieldNames),
                    SAPFeats.(SAPFeat_FieldNames{k}) = [SAPFeats.(SAPFeat_FieldNames{k}); Feats.(SAPFeat_FieldNames{k})(:)];
                end
            end
        else
            Song = [];
            for Files = BoutOnsetFile:BoutOffsetFile,
                [Temp, Fs] = GetData(BirdParameters.DataDirectory, BirdParameters.SongFileNames{Files}, BirdParameters.FileType, 0);
                Song = [Song; Temp(:)];
            end
            
            Time = (1:1:length(Song))/Fs;
            
            Onsets = BirdParameters.ActualAllOnsets(BirdParameters.BoutIndices(j,1):BirdParameters.BoutIndices(j,2)) - ContinuousFileTime(BoutOnsetFile);
            Offsets = BirdParameters.ActualAllOffsets(BirdParameters.BoutIndices(j,1):BirdParameters.BoutIndices(j,2)) - ContinuousFileTime(BoutOnsetFile);
                
            SongBoutOnsetIndex = find(Time >= (Onsets(1)/1000 - Padding), 1, 'first');
            SongBoutOffsetIndex = find(Time >= (Offsets(end)/1000 + Padding), 1, 'first');
                
            [Feats, RawFeats, FeatsFs] = ASSLCalculateSAPFeatsWithOnsets(Song(SongBoutOnsetIndex:SongBoutOffsetIndex), Time(SongBoutOnsetIndex:SongBoutOffsetIndex), Fs, Onsets/1000, Offsets/1000);
            if (~exist('SAPFeat_FieldNames', 'var'))
                SAPFeat_FieldNames = fieldnames(Feats);
                for k = 1:length(SAPFeat_FieldNames),
                    SAPFeats.(SAPFeat_FieldNames{k}) = Feats.(SAPFeat_FieldNames{k})(:);
                end
            else
                for k = 1:length(SAPFeat_FieldNames),
                    SAPFeats.(SAPFeat_FieldNames{k}) = [SAPFeats.(SAPFeat_FieldNames{k}); Feats.(SAPFeat_FieldNames{k})(:)];
                end
            end
        end
    end
end
fprintf('\n');