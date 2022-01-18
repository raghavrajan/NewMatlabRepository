function [BirdParameters, Flag] = ProcessSongData(BirdDetailsTextFile, InterboutInterval)

% First get details from the CSV text file
disp('Getting header data from CSV file ...');
[HeaderLine, BirdDetails] = LSINA_GetDetailsFromCSVFile(BirdDetailsTextFile);

% Now parse all the lines into the appropriate variables based on the
% header line
disp('Getting data from CSV file ...');
[BirdParameters] = LSINA_ParseHeaderBirdData(HeaderLine, BirdDetails);

% Now put in the input inter-bout interval
for i = 1:length(BirdParameters),
    BirdParameters(i).Interboutinterval = InterboutInterval;
end

% Now for each of the birds, load up all the filenames
disp('Loading up filenames ...');
for i = 1:length(BirdParameters),
    FileNames_FileName = [BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.', BirdParameters(i).Condition, '.', BirdParameters(i).Microphone, '.FileNames.mat'];
    if (exist(FileNames_FileName, 'file'))
        load(FileNames_FileName);
        BirdParameters(i).SongFileNames = SongFileNames;
    else
        [BirdParameters(i).SongFileNames] = LSINA_GetDataFileNames(BirdParameters(i));
        SongFileNames = BirdParameters(i).SongFileNames;
        save(FileNames_FileName, 'SongFileNames');
    end
end

% Now check if notes files have unequal numbers of labels, onsets and
% offsets
disp('Checking note data ...');
for i = 1:length(BirdParameters),
    if (exist(BirdParameters(i).DataDirectory, 'dir'))
        if (isfield(BirdParameters(i), 'UndirSongFileList'))
            Flag(i) = CheckLengthsOnsetsOffsetsLabels(BirdParameters(i).DataDirectory, BirdParameters(i).UndirSongFileList);
        else
            Flag(i) = CheckLengthsOnsetsOffsetsLabels(BirdParameters(i).DataDirectory, BirdParameters(i).SongFileList);
        end
    else
        Flag(i) = 0;
    end
end

% If there are files with unequal numbers, then don't run the script - just
% return
if (sum(Flag) > 0)
    return;
else
    disp('All notes files are ok');
end


% Now load up the note files and the length of each file
disp('Loading up note data ...');
for i = 1:length(BirdParameters),
    Notes_FileName = [BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.', BirdParameters(i).Condition, '.', BirdParameters(i).Microphone, '.Notes.mat'];
    if (exist(Notes_FileName, 'file'))
        load(Notes_FileName);
        BirdParameters(i).NoteInfo = NoteInfo;
        BirdParameters(i).FileLen = FileLen;
    else
        [BirdParameters(i).NoteInfo, BirdParameters(i).FileLen] = LSINA_LoadNoteFileInfo(BirdParameters(i));
        NoteInfo = BirdParameters(i).NoteInfo;
        FileLen = BirdParameters(i).FileLen;
        save(Notes_FileName, 'NoteInfo', 'FileLen');
    end
end

% % Now calculate template match values for each of the syllables
% disp('Calculating template match values for all syllables ...');
% for i = 1:length(BirdParameters),
%     TemplateMatchValues_FileName = [BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.', BirdParameters(i).Condition, '.', BirdParameters(i).Microphone, '.TemplateMatchValues.mat'];
%     if (exist(TemplateMatchValues_FileName, 'file'))
%         load(TemplateMatchValues_FileName);
%         BirdParameters(i).TemplateMatchValues = TemplateMatchValues;
%     else
%         [BirdParameters(i).TemplateMatchValues] = LSINA_CalcTemplateMatchValues(BirdParameters(i));
%         TemplateMatchValues = BirdParameters(i).TemplateMatchValues;
%         save(TemplateMatchValues_FileName, 'TemplateMatchValues');
%     end
% end

% % Now calculate log amplitudes based on Aronov Fee and take into
% % consideration the baseline.
% disp('Calculating log amplitude, Aronov, Fee for all syllables ...');
% for i = 1:length(BirdParameters),
%     AronovFeeLogAmpValues_FileName = [BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.', BirdParameters(i).Condition, '.', BirdParameters(i).Microphone, '.AronovFeeLogAmpValues.mat'];
%     if (exist(AronovFeeLogAmpValues_FileName, 'file'))
%         load(AronovFeeLogAmpValues_FileName);
%         BirdParameters(i).AronovFeeLogAmpValues = AronovFeeLogAmpValues;
%         BirdParameters(i).BaselineAmpValue = BaselineAmpValue;
%         BirdParameters(i).SyllOnsetValues = SyllOnsetValues;
%         BirdParameters(i).SyllOffsetValues = SyllOffsetValues;
%         BirdParameters(i).FirstHalfAmpValue = FirstHalfAmpValue;
%         BirdParameters(i).SecondHalfAmpValue = SecondHalfAmpValue;
%     else
%         [BirdParameters(i).AronovFeeLogAmpValues, BirdParameters(i).BaselineAmpValue, BirdParameters(i).SyllOnsetValues, BirdParameters(i).SyllOffsetValues, BirdParameters(i).FirstHalfAmpValue, BirdParameters(i).SecondHalfAmpValue] = LSINA_CalcAronovFeeLogAmpValues(BirdParameters(i));
%         AronovFeeLogAmpValues = BirdParameters(i).AronovFeeLogAmpValues;
%         BaselineAmpValue = BirdParameters(i).BaselineAmpValue;
%         SyllOnsetValues = BirdParameters(i).SyllOnsetValues;
%         SyllOffsetValues = BirdParameters(i).SyllOffsetValues;
%         FirstHalfAmpValue = BirdParameters(i).FirstHalfAmpValue;
%         SecondHalfAmpValue = BirdParameters(i).SecondHalfAmpValue;
%         
%         save(AronovFeeLogAmpValues_FileName, 'AronovFeeLogAmpValues', 'BaselineAmpValue', 'SyllOnsetValues', 'SyllOffsetValues', 'FirstHalfAmpValue', 'SecondHalfAmpValue');
%     end
% end

% Now calculate log amplitudes based on SAP (deriv.m) and also note the
% baseline so that I can subtract that later
% For the moment am going to allow SAP code to subtract -70dB as per its
% default
disp('Calculating log amplitude, SAP for all syllables ...');
for i = 1:length(BirdParameters),
    if (mod(i,5) == 0)
        fprintf('>%i>', i);
    end
    FFTLogAmpValues_FileName = [BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.', BirdParameters(i).Condition, '.', BirdParameters(i).Microphone, '.FFTLogAmpValues.mat'];
    if (exist(FFTLogAmpValues_FileName, 'file'))
        load(FFTLogAmpValues_FileName);
        BirdParameters(i).FFTLogAmplitudes = FFTLogAmplitudes;
        BirdParameters(i).BaselineAmpValue = BaselineAmpValue;
        BirdParameters(i).SyllableStatus = SyllableStatus;
        BirdParameters(i).Amp_Fs = Amp_Fs;
    else
        [BirdParameters(i).FFTLogAmplitudes, BirdParameters(i).BaselineAmpValue, BirdParameters(i).SyllableStatus, BirdParameters(i).Amp_Fs] = LSINA_CalcFFTLogAmpValues(BirdParameters(i));
        FFTLogAmplitudes = BirdParameters(i).FFTLogAmplitudes;
        BaselineAmpValue = BirdParameters(i).BaselineAmpValue;
        SyllableStatus = BirdParameters(i).SyllableStatus;
        Amp_Fs = BirdParameters(i).Amp_Fs;
        save(FFTLogAmpValues_FileName, 'FFTLogAmplitudes', 'BaselineAmpValue', 'SyllableStatus', 'Amp_Fs');
    end
end
fprintf('\n');

% Now the first thing to do would be to check if it is continuous data or
% not. If it is continuous data, then check for consecutive Capital letter
% syllables that would correspond to the same syllable split over two
% conseecutive files. These have to be merged.
% Do this for continuous data and then put together one long list of
% syllables and their corresponding file #s, onsets and offsets.

disp('Putting together list of syllables ...');
for i = 1:length(BirdParameters),
    SyllableInfo_FileName = [BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.', BirdParameters(i).Condition, '.', BirdParameters(i).Microphone, '.SyllableInfo.mat'];
    if (exist(SyllableInfo_FileName, 'file'))
        load(SyllableInfo_FileName);
        BirdParameters(i).SyllableData = SyllableData;
    else
        [BirdParameters(i).SyllableData] = GetSyllableListInfo(BirdParameters(i));
        SyllableData = BirdParameters(i).SyllableData;
        save(SyllableInfo_FileName, 'SyllableData');
    end
    fprintf('%d >> ', i);
end
fprintf('\n');

% Now calculate SAP features for all the files
disp('Calculating SAP features ...');
for i = 1:length(BirdParameters),
    disp(['Calculating SAP features for ', BirdParameters(i).BirdName, '-', BirdParameters(i).DataLabel, ' ...']);
    SAPFileName = [BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.', BirdParameters(i).Condition, '.', BirdParameters(i).Microphone, '.SAPFeats.mat'];
    if (exist(SAPFileName, 'file'))
        load(SAPFileName);
        BirdParameters(i).SAPFeatsMatrix = SAPFeatsMatrix;
        BirdParameters(i).SAPFeat_FieldNames = SAPFeat_FieldNames;
    else
        [BirdParameters(i).SAPFeatsMatrix, BirdParameters(i).SAPFeat_FieldNames] = CalcSAPFeats(BirdParameters(i));
        [FFValues] = CalcFF(BirdParameters(i));
        BirdParameters(i).SAPFeatsMatrix(:,end) = FFValues(:,2);
        BirdParameters(i).SAPFeatsMatrix(:,end+1) = FFValues(:,1);
        BirdParameters(i).SAPFeat_FieldNames{end} = 'MeanFundamentalFrequency';
        BirdParameters(i).SAPFeat_FieldNames{end+1} = 'MedianFundamentalFrequency';
        SAPFeatsMatrix = BirdParameters(i).SAPFeatsMatrix;
        SAPFeat_FieldNames = BirdParameters(i).SAPFeat_FieldNames;
        save(SAPFileName, 'SAPFeatsMatrix', 'SAPFeat_FieldNames');
    end        
    % Now, I need to remove outliers for each of the syllables.
    % For each syllable type I will calculate upper and lower outlier
    % threshold as (75th percentile + 3*inter-quartile range) and (25th
    % percentile - 3*inter-quartile range) respectively. For each
    % syllable type, then I will check if any of the syllables are
    % above or below the upper and lower threshold. If they are then I
    % will make them NaN

%     AmplitudeColIndex = find(cellfun(@length, strfind(BirdParameters(i).SAPFeat_FieldNames, 'LogAmplitude')));
%     UniqueSylls = unique(char(BirdParameters(i).SyllableData(:,1)));
%     for j = 1:length(UniqueSylls),
%         Indices = find(char(BirdParameters(i).SyllableData(:,1)) == UniqueSylls(j));
%         if (~isempty(Indices))
%             % for k = 1:size(BirdParameters(i).SAPFeatsMatrix,2),
%             for k = AmplitudeColIndex, % removing outliers only for amplitude
%                 UpperOutlierThreshold = prctile(BirdParameters(i).SAPFeatsMatrix(Indices,k), 75) + 3*iqr(BirdParameters(i).SAPFeatsMatrix(Indices,k));
%                 LowerOutlierThreshold = prctile(BirdParameters(i).SAPFeatsMatrix(Indices,k), 25) - 3*iqr(BirdParameters(i).SAPFeatsMatrix(Indices,k));
% 
%                 OutlierIndices = find((BirdParameters(i).SAPFeatsMatrix(Indices,k) > UpperOutlierThreshold) | (BirdParameters(i).SAPFeatsMatrix(Indices,k) < LowerOutlierThreshold));
%                 disp([BirdParameters(i).BirdName, ': Syllable ', UniqueSylls(j), '; removed ', num2str(length(OutlierIndices)), ' outliers out of a total of ', num2str(length(Indices)), '(', num2str(100*length(OutlierIndices)/length(Indices)), '%)']);
%                 % If a syllable is an outlier in any one of the
%                 % features, I make all its feature values as NaN,
%                 % effectively excluding that syllable completely.
%                 BirdParameters(i).SAPFeatsMatrix(Indices(OutlierIndices), :) = NaN;
%             end
%         end
%     end
end

% Now calculate rms amplitue
disp('Calculating rms amplitude ...');
for i = 1:length(BirdParameters),
    disp(['Calculating rms amplitude  for ', BirdParameters(i).BirdName, '-', BirdParameters(i).DataLabel, ' ...']);
    SAPFileName = [BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.', BirdParameters(i).Condition, '.', BirdParameters(i).Microphone, '.AmplitudeRMS.mat'];
    if (exist(SAPFileName, 'file'))
        load(SAPFileName);
        BirdParameters(i).AmplitudeRMS = AmplitudeRMS;
    else
        [BirdParameters(i).AmplitudeRMS] = CalcFeats_AmplitudeRMS(BirdParameters(i));
        AmplitudeRMS = BirdParameters(i).AmplitudeRMS;
        save(SAPFileName, 'AmplitudeRMS');
    end 
end

% Now calculate log amplitude Kao et al. 2005
disp('Calculating amplitude Kao ...');
for i = 1:length(BirdParameters),
    disp(['Calculating amplitude Kao for ', BirdParameters(i).BirdName, '-', BirdParameters(i).DataLabel, ' ...']);
    SAPFileName = [BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.', BirdParameters(i).Condition, '.', BirdParameters(i).Microphone, '.AmplitudeKao.mat'];
    if (exist(SAPFileName, 'file'))
        load(SAPFileName);
        BirdParameters(i).AmplitudeKao = AmplitudeKao;
    else
        [BirdParameters(i).AmplitudeKao] = CalcFeats_AmplitudeKao(BirdParameters(i));
        AmplitudeKao = BirdParameters(i).AmplitudeKao;
        save(SAPFileName, 'AmplitudeKao');
    end        
%     % Now, I need to remove outliers for each of the syllables.
%     % For each syllable type I will calculate upper and lower outlier
%     % threshold as (75th percentile + 3*inter-quartile range) and (25th
%     % percentile - 3*inter-quartile range) respectively. For each
%     % syllable type, then I will check if any of the syllables are
%     % above or below the upper and lower threshold. If they are then I
%     % will make them NaN
% 
%     UniqueSylls = unique(char(BirdParameters(i).SyllableData(:,1)));
%     for j = 1:length(UniqueSylls),
%         Indices = find(char(BirdParameters(i).SyllableData(:,1)) == UniqueSylls(j));
%         if (~isempty(Indices))
%             UpperOutlierThreshold = prctile(BirdParameters(i).AmplitudeKao(Indices), 75) + 3*iqr(BirdParameters(i).AmplitudeKao(Indices));
%             LowerOutlierThreshold = prctile(BirdParameters(i).AmplitudeKao(Indices), 25) - 3*iqr(BirdParameters(i).AmplitudeKao(Indices));
% 
%             OutlierIndices = find((BirdParameters(i).AmplitudeKao(Indices) > UpperOutlierThreshold) | (BirdParameters(i).AmplitudeKao(Indices) < LowerOutlierThreshold));
%             disp([BirdParameters(i).BirdName, ': Syllable ', UniqueSylls(j), '; removed ', num2str(length(OutlierIndices)), ' outliers out of a total of ', num2str(length(Indices)), '(', num2str(100*length(OutlierIndices)/length(Indices)), '%), amplitude Kao']);
%             % If a syllable is an outlier in any one of the
%             % features, I make all its feature values as NaN,
%             % effectively excluding that syllable completely.
%             BirdParameters(i).AmplitudeKao(Indices(OutlierIndices)) = NaN;
%         end
%     end
end

% Now calculate FF Kao et al. 2005 using autocorr
disp('Calculating FF autocorr Kao ...');
for i = 1:length(BirdParameters),
    disp(['Calculating FF autocorr Kao for ', BirdParameters(i).BirdName, '-', BirdParameters(i).DataLabel, ' ...']);
    SAPFileName = [BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.', BirdParameters(i).Condition, '.', BirdParameters(i).Microphone, '.FFAutocorrKao.mat'];
    if (exist(SAPFileName, 'file'))
        load(SAPFileName);
        BirdParameters(i).FFAutocorrKao = FFAutocorrKao;
    else
        [BirdParameters(i).FFAutocorrKao] = CalcFF_AutoCorr(BirdParameters(i));
        % BirdParameters(i).FFAutocorrKao = ones(size(BirdParameters(i).AmplitudeKao))*NaN;
        FFAutocorrKao = BirdParameters(i).FFAutocorrKao;
        save(SAPFileName, 'FFAutocorrKao');
    end        
end

% Now split up the files into bouts based on inter-bout interval that is
% also specified in the .csv file
disp('Identifying bouts ...');
for i = 1:length(BirdParameters),
    fprintf('%d >> ', i);
    BoutInfo_FileName = [BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.', BirdParameters(i).Condition, '.', BirdParameters(i).Microphone, '.', num2str(BirdParameters(i).Interboutinterval), '.BoutInfo.mat'];
    if (exist(BoutInfo_FileName, 'file'))
        load(BoutInfo_FileName);
        BirdParameters(i).Bouts = Bouts;
    else
        [BirdParameters(i).Bouts] = LSINA_GetBoutInfo(BirdParameters(i));
        Bouts = BirdParameters(i).Bouts;
        save(BoutInfo_FileName, 'Bouts');
    end
    fprintf('%d >> ', i);
end
fprintf('\n');

% Now to write the bout information into a text file - I want to write bout
% #, filename, bout onset and offset in sec to a text file
disp('Writing bout info to text file ...');
for i = 1:length(BirdParameters),
    fprintf('%d >> ', i);
    Fid = fopen([BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.', BirdParameters(i).Condition, '.', BirdParameters(i).Microphone, '.BoutInfo.txt'], 'w');
    fprintf(Fid, 'Bout #\tFileName\tBout onset (sec)\tBout offset (sec)\n');
    for j = 1:size(BirdParameters(i).Bouts, 1),
        if ((BirdParameters(i).Bouts(j,7) == 1) && (BirdParameters(i).Bouts(j,8) > 0) && (BirdParameters(i).Bouts(j,9) > 1))
            fprintf(Fid, '%d\t%s\t%3.5f\t%3.5f\n', j, BirdParameters(i).SongFileNames{BirdParameters(i).Bouts(j,3)}, BirdParameters(i).Bouts(j,5)/1000, BirdParameters(i).Bouts(j,6)/1000);
        end
    end
    fclose(Fid);
end
fprintf('\n');

disp('Finished Analysis');