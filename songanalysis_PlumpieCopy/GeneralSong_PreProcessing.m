function [SURecordingDetails] = GeneralSong_PreProcessing(SongDetailsFile, InterBoutInterval, SavedFileDir)

% ======== Housekeeping stuff common to all figures =======================

% Some common parameters

% First, read in the csv file and get all the columns into appropriate
% fields of a structure called RecordingDataParameters

% First for single units
if (exist(SongDetailsFile, 'file'))
    disp('Getting header data from CSV file ...');
    [HeaderLine, BirdDetails] = LSINA_GetDetailsFromCSVFile(SongDetailsFile);

    disp('Getting data from CSV file ...');
    [SURecordingDetails] = LSINA_ParseHeaderBirdData(HeaderLine, BirdDetails);
end

if (exist('SURecordingDetails', 'var'))
    % If saved file exists, then load that up
    for i = 1:length(SURecordingDetails),
        disp(['Song file list #', num2str(i), ' ...']);
        
        %SURecordingDetails(i).Interboutinterval = InterBoutInterval;
        
        % First load up all the things that are independent of the
        % inter-bout interval
        
        disp('Loading up filenames ...');
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.FileNames.mat']);
        if (exist(SavedFileName, 'file'))
            load(SavedFileName);
            SURecordingDetails(i).SongFileNames = Temp(1).SongFileNames;
        else
            [SURecordingDetails(i).SongFileNames] = LSINA_GetDataFileNames(SURecordingDetails(i));
            clear Temp;
            Temp(1).SongFileNames = SURecordingDetails(i).SongFileNames;
            save(SavedFileName, 'Temp');
        end
        
        % Now get the labels, spike data from the respective files for each of
        % the files
        disp('Loading up note data ...');
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.NoteInfo.mat']);
        if (exist(SavedFileName, 'file'))
            load(SavedFileName);
            SURecordingDetails(i).ActualNoteInfo = Temp(1).ActualNoteInfo;
            SURecordingDetails(i).NoteInfo = Temp(1).NoteInfo;
            SURecordingDetails(i).FileLen = Temp(1).FileLen;
        else
            [SURecordingDetails(i).NoteInfo, SURecordingDetails(i).FileLen] = LSINA_LoadNoteFileInfo(SURecordingDetails(i));
            
            % Now copy noteinfo to another field called ActualNoteInfo and
            % remove the female in, out, door close event labels and female
            % call labels from NoteInfo
            SURecordingDetails(i).ActualNoteInfo = SURecordingDetails(i).NoteInfo;
            for j = 1:length(SURecordingDetails(i).NoteInfo),
                if (~isempty(SURecordingDetails(i).NoteInfo{j}.labels))
                    FemaleINOUTEvents = regexp(char(SURecordingDetails(i).NoteInfo{j}.labels), ['[', SURecordingDetails(i).Femalecalllabel, SURecordingDetails(i).FemaleinLabel, SURecordingDetails(i).FemaleoutLabel, SURecordingDetails(i).Closedoorlabel, ']']);
                    SURecordingDetails(i).NoteInfo{j}.labels(FemaleINOUTEvents) = [];
                    SURecordingDetails(i).NoteInfo{j}.onsets(FemaleINOUTEvents) = [];
                    SURecordingDetails(i).NoteInfo{j}.offsets(FemaleINOUTEvents) = [];
                end
            end
        
            clear Temp;
            Temp(1).NoteInfo = SURecordingDetails(i).NoteInfo;
            Temp(1).ActualNoteInfo = SURecordingDetails(i).ActualNoteInfo;
            Temp(1).FileLen = SURecordingDetails(i).FileLen;
            save(SavedFileName, 'Temp');
        end
        
        % Now the first thing to do would be to check if it is continuous data or
        % not. If it is continuous data, then check for consecutive Capital letter
        % syllables that would correspond to the same syllable split over two
        % conseecutive files. These have to be merged.
        % Do this for continuous data and then put together one long list of
        % syllables and their corresponding file #s, onsets and offsets.

        disp('Putting together list of syllables ...');
        Saved_FileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.SyllableInfo.mat']);
        if (exist(Saved_FileName, 'file'))
            load(Saved_FileName);
            SURecordingDetails(i).SyllableData = SyllableData;
        else
            [SURecordingDetails(i).SyllableData] = GetSyllableListInfo(SURecordingDetails(i));
            SyllableData = SURecordingDetails(i).SyllableData;
            save(Saved_FileName, 'SyllableData');
        end

        % Now calculate SAP features for all the files
        disp('Calculating SAP features ...');
        disp(['Calculating SAP features for ', SURecordingDetails(i).BirdName, '-', SURecordingDetails(i).DataLabel, ' ...']);
        Saved_FileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.SAPFeats.mat']);
        if (exist(Saved_FileName, 'file'))
            load(Saved_FileName);
            SURecordingDetails(i).SAPFeatsMatrix = SAPFeatsMatrix;
            SURecordingDetails(i).SAPFeat_FieldNames = SAPFeat_FieldNames;
        else
            [SURecordingDetails(i).SAPFeatsMatrix, SURecordingDetails(i).SAPFeat_FieldNames] = CalcSAPFeats(SURecordingDetails(i));
            [FFValues] = CalcFF(SURecordingDetails(i));
            SURecordingDetails(i).SAPFeatsMatrix(:,end) = FFValues(:,2);
            SURecordingDetails(i).SAPFeatsMatrix(:,end+1) = FFValues(:,1);
            SURecordingDetails(i).SAPFeat_FieldNames{end} = 'MeanFundamentalFrequency';
            SURecordingDetails(i).SAPFeat_FieldNames{end+1} = 'MedianFundamentalFrequency';

%             % Now, I need to remove outliers for each of the syllables.
%             % For each syllable type I will calculate upper and lower outlier
%             % threshold as (75th percentile + 3*inter-quartile range) and (25th
%             % percentile - 3*inter-quartile range) respectively. For each
%             % syllable type, then I will check if any of the syllables are
%             % above or below the upper and lower threshold. If they are then I
%             % will make them NaN
%             UniqueSylls = unique(char(SURecordingDetails(i).SyllableData(:,1)));
%             for j = 1:length(UniqueSylls),
%                 Indices = find(char(SURecordingDetails(i).SyllableData(:,1)) == UniqueSylls(j));
%                 if (~isempty(Indices))
%                     for k = 1:size(SURecordingDetails(i).SAPFeatsMatrix,2),
%                         UpperOutlierThreshold = prctile(SURecordingDetails(i).SAPFeatsMatrix(Indices,k), 75) + 3*iqr(SURecordingDetails(i).SAPFeatsMatrix(Indices,k));
%                         LowerOutlierThreshold = prctile(SURecordingDetails(i).SAPFeatsMatrix(Indices,k), 25) - 3*iqr(SURecordingDetails(i).SAPFeatsMatrix(Indices,k));
% 
%                         OutlierIndices = find((SURecordingDetails(i).SAPFeatsMatrix(Indices,k) > UpperOutlierThreshold) | (SURecordingDetails(i).SAPFeatsMatrix(Indices,k) < LowerOutlierThreshold));
%                         disp([SURecordingDetails(i).BirdName, ': Syllable ', UniqueSylls(j), '; removed ', num2str(length(OutlierIndices)), ' outliers out of a total of ', num2str(length(Indices)), '(', num2str(100*length(OutlierIndices)/length(Indices)), '%)']);
%                         % If a syllable is an outlier in any one of the
%                         % features, I make all its feature values as NaN,
%                         % effectively excluding that syllable completely.
%                         SURecordingDetails(i).SAPFeatsMatrix(Indices(OutlierIndices), :) = NaN;
%                     end
%                 end
%             end

            SAPFeatsMatrix = SURecordingDetails(i).SAPFeatsMatrix;
            SAPFeat_FieldNames = SURecordingDetails(i).SAPFeat_FieldNames;
            save(Saved_FileName, 'SAPFeatsMatrix', 'SAPFeat_FieldNames');
        end
        
        % Now split up the files into bouts based on inter-bout interval that is
        % also specified in the .csv file
        disp('Identifying bouts ...');
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.InterboutInterval', num2str(InterBoutInterval), '.BoutInfo.mat']);
        if (exist(SavedFileName, 'file'))
            load(SavedFileName);
            SURecordingDetails(i).Bouts = Temp(1).Bouts;
            SURecordingDetails(i).Gaps = Temp(1).Gaps;
        else
            [SURecordingDetails(i).Bouts, SURecordingDetails(i).Gaps] = LSINA_GetBoutInfo(SURecordingDetails(i));
            clear Temp;
            Temp(1).Bouts = SURecordingDetails(i).Bouts;
            Temp(1).Gaps = SURecordingDetails(i).Gaps;
            save(SavedFileName, 'Temp');
        end
        
        % Now find out where the directed presentations are - starting file and
        % starting time, ending file and ending time
        disp('Locating directed presentations ...');
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.InterboutInterval', num2str(InterBoutInterval), '.DirPresentationInfo.mat']);
        if (exist(SavedFileName, 'file'))
            load(SavedFileName);
            SURecordingDetails(i).DirPresentations = Temp(1).DirPresentations;
        else
            [SURecordingDetails(i).DirPresentations] = LSINA_LocateDirPresentations(SURecordingDetails(i));
            clear Temp;
            Temp(1).DirPresentations = SURecordingDetails(i).DirPresentations;
            save(SavedFileName, 'Temp');
        end
            
        % Now check for each bout whether it is dir or undir
        disp('Labelling bouts as dir (1) or undir (0) ...');
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.InterboutInterval', num2str(InterBoutInterval), '.BoutDirUndirInfo.mat']);
        if (exist(SavedFileName, 'file'))
            load(SavedFileName);
            SURecordingDetails(i).BoutDirUnDir = Temp(1).BoutDirUnDir;
        else
            [SURecordingDetails(i).BoutDirUnDir] = LSINA_LabelBoutsDirUnDir(SURecordingDetails(i));
            clear Temp;
            Temp(1).BoutDirUnDir = SURecordingDetails(i).BoutDirUnDir;
            save(SavedFileName, 'Temp');
        end
            
        % Now check for each gap whether it is dir or undir
        disp('Labelling gaps as dir (1) or undir (0) ...');
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.InterboutInterval', num2str(InterBoutInterval), '.GapDirUndirInfo.mat']);
        if (exist(SavedFileName, 'file'))
            load(SavedFileName);
            SURecordingDetails(i).GapDirUnDir = Temp(1).GapDirUnDir;
        else
            [SURecordingDetails(i).Gaps, SURecordingDetails(i).GapDirUnDir] = LSINA_LabelGapsDirUnDir(SURecordingDetails(i));
            clear Temp;
            Temp(1).GapDirUnDir = SURecordingDetails(i).GapDirUnDir;
            save(SavedFileName, 'Temp');
        end
            
        % Now for each neuron, find all the bouts that have the required
        % inter-bout interval pre and post and then get the spike trains
        % Also, get the spontaneous window spike count, the pre-bout window
        % spike count, the PST and the raster for each trial. Also, get the
        % identities of the first syllable. 
        % This can all be saved to a .mat file for faster loading and for
        % all further analysis
        disp('Finding dir and undir bouts ...');
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.InterboutInterval', num2str(InterBoutInterval), '.BoutDetails.mat']);
        if (exist(SavedFileName, 'file'))
            load(SavedFileName);
            SURecordingDetails(i).UnDirBoutDetails = Temp(1).UnDirBoutDetails;
            SURecordingDetails(i).DirBoutDetails = Temp(1).DirBoutDetails;
        else
            [SURecordingDetails(i).UnDirBoutDetails, SURecordingDetails(i).DirBoutDetails] = LSINA_FindDirUnDirBouts(SURecordingDetails(i));        
            clear Temp;
            Temp(1).UnDirBoutDetails = SURecordingDetails(i).UnDirBoutDetails;
            Temp(1).DirBoutDetails = SURecordingDetails(i).DirBoutDetails;
            save(SavedFileName, 'Temp');
        end
        
        % Now for each neuron, make a plot which shows each file as one
        % row and shows the onsets and offsets of bouts coloured as
        % blue or red or black (depending on whether it is undir, dir
        % or ambiguous respectively) as the next row, and then same for
        % gaps as the next row and then finally the last row showing
        % the onset and offset of dir presentations. THis will help me
        % verify that the bouts and gaps have been correctly labelled
        % as dir/undir/ambiguous.
        %disp('Showing bout and gap structure in each file ...');
        %LSINA_ShowBoutGapStructureForEachFile(SURecordingDetails(i));
        
        % Now for each neuron, go through all the note files and see if the
        % lengths of onsets, offsets and labels are not equal for any of
        % the files
        disp('Checking note files');
        for j = 1:length(SURecordingDetails(i).NoteInfo),
            if ((length(SURecordingDetails(i).NoteInfo{j}.onsets) ~= length(SURecordingDetails(i).NoteInfo{j}.offsets)) || (length(SURecordingDetails(i).NoteInfo{j}.labels) ~= length(SURecordingDetails(i).NoteInfo{j}.offsets)) || (length(SURecordingDetails(i).NoteInfo{j}.onsets) ~= length(SURecordingDetails(i).NoteInfo{j}.labels)))
                disp(SURecordingDetails(i).SongFileNames{j});
            end
        end
    end
end
    
