function [SURecordingDetails] = Common_NeuralAnalysis_PreProcessing(InterBoutInterval, SingleUnitDataFile, MultiUnitDataFile, SavedFileDir)

% ======== Housekeeping stuff common to all figures =======================

% First, read in the csv file and get all the columns into appropriate
% fields of a structure called RecordingDataParameters

% First for single units
if (exist(SingleUnitDataFile, 'file'))
    disp('Getting header data from CSV file ...');
    [HeaderLine, BirdDetails] = LSINA_GetDetailsFromCSVFile(SingleUnitDataFile);

    disp('Getting data from CSV file ...');
    [SURecordingDetails] = LSINA_ParseHeaderBirdData(HeaderLine, BirdDetails);
end

% Then for multi units
if (exist(MultiUnitDataFile, 'file'))
    MURecordingDataParameters = ParseRecordingDataFile(MultiUnitDataFile);
end

if (exist('SURecordingDetails', 'var'))
    % If saved file exists, then load that up
    for i = 1:length(SURecordingDetails),
        disp(['Neuron #', num2str(i), ' ...']);
        
        SURecordingDetails(i).Interboutinterval = InterBoutInterval;
        
        % First load up all the things that are independent of the
        % inter-bout interval
        
        disp('Loading up filenames ...');
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.', SURecordingDetails(i).Neurontype, '.SpikeChanNo', num2str(SURecordingDetails(i).SpikeChanNo), '.FileNames.mat']);
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
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.', SURecordingDetails(i).Neurontype, '.SpikeChanNo', num2str(SURecordingDetails(i).SpikeChanNo), '.NoteInfo.mat']);
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
        
        
        disp('Loading up spike data ...');
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.', SURecordingDetails(i).Neurontype, '.SpikeChanNo', num2str(SURecordingDetails(i).SpikeChanNo), '.SpikeInfo.mat']);
        if (exist(SavedFileName, 'file'))
            load(SavedFileName);
            SURecordingDetails(i).SpikeInfo = Temp(1).SpikeInfo;
        else
            [SURecordingDetails(i).SpikeInfo] = LSINA_LoadSpikeTimes(SURecordingDetails(i));
            clear Temp;
            Temp(1).SpikeInfo = SURecordingDetails(i).SpikeInfo;
            save(SavedFileName, 'Temp');
        end
        
        disp('Now getting spikes from different clusters ...');
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.', SURecordingDetails(i).Neurontype, '.SpikeChanNo', num2str(SURecordingDetails(i).SpikeChanNo), '.SpikeClusterInfo.mat']);
        if (exist(SavedFileName, 'file'))
            load(SavedFileName);
            SURecordingDetails(i).ClusterSpikeTimes = Temp(1).ClusterSpikeTimes;
            SURecordingDetails(i).OptionalClusterSpikeTimes = Temp(1).OptionalClusterSpikeTimes;
        else
            [SURecordingDetails(i).ClusterSpikeTimes, SURecordingDetails(i).OptionalClusterSpikeTimes] = LSINA_GetClusterSpikeTimes(SURecordingDetails(i));
            clear Temp;
            Temp(1).ClusterSpikeTimes = SURecordingDetails(i).ClusterSpikeTimes;
            Temp(1).OptionalClusterSpikeTimes = SURecordingDetails(i).OptionalClusterSpikeTimes;
            save(SavedFileName, 'Temp');
        end
        
        % Now split up the files into bouts based on inter-bout interval that is
        % also specified in the .csv file
        disp('Identifying bouts ...');
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.', SURecordingDetails(i).Neurontype, '.InterboutInterval', num2str(InterBoutInterval), '.SpikeChanNo', num2str(SURecordingDetails(i).SpikeChanNo), '.BoutInfo.mat']);
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
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.', SURecordingDetails(i).Neurontype, '.InterboutInterval', num2str(InterBoutInterval), '.SpikeChanNo', num2str(SURecordingDetails(i).SpikeChanNo), '.DirPresentationInfo.mat']);
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
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.', SURecordingDetails(i).Neurontype, '.InterboutInterval', num2str(InterBoutInterval), '.SpikeChanNo', num2str(SURecordingDetails(i).SpikeChanNo), '.BoutDirUndirInfo.mat']);
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
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.', SURecordingDetails(i).Neurontype, '.InterboutInterval', num2str(InterBoutInterval), '.SpikeChanNo', num2str(SURecordingDetails(i).SpikeChanNo), '.GapDirUndirInfo.mat']);
        if (exist(SavedFileName, 'file'))
            load(SavedFileName);
            SURecordingDetails(i).GapDirUnDir = Temp(1).GapDirUnDir;
        else
            [SURecordingDetails(i).Gaps, SURecordingDetails(i).GapDirUnDir] = LSINA_LabelGapsDirUnDir(SURecordingDetails(i));
            clear Temp;
            Temp(1).GapDirUnDir = SURecordingDetails(i).GapDirUnDir;
            save(SavedFileName, 'Temp');
        end
            
        % Now get the waveforms for all the spikes for each neuron,
        % align them and then get the SNR for that neuron.
        disp('Getting spike waveforms and SNR for each neuron ...');
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.', SURecordingDetails(i).Neurontype, '.InterboutInterval', num2str(InterBoutInterval), '.SpikeChanNo', num2str(SURecordingDetails(i).SpikeChanNo), '.SpikeWaveformsSNR.mat']);
        if (exist(SavedFileName, 'file'))
            load(SavedFileName);
            SURecordingDetails(i).SpikeWaveforms = Temp(1).SpikeWaveforms;
            % SURecordingDetails(i).AlignedSpikeWaveforms = Temp(1).AlignedSpikeWaveforms;
            SURecordingDetails(i).SNR = Temp(1).SNR;
        else
            [SURecordingDetails(i).SpikeWaveforms, SURecordingDetails(i).SNR] = LSINA_GetSpikeWaveforms_SNR(SURecordingDetails(i));            
            clear Temp;
            Temp(1).SpikeWaveforms = SURecordingDetails(i).SpikeWaveforms;
            % Temp(1).AlignedSpikeWaveforms = SURecordingDetails(i).AlignedSpikeWaveforms;
            Temp(1).SNR = SURecordingDetails(i).SNR;
            save(SavedFileName, 'Temp');
        end

        % Now for each neuron, find all the bouts that have the required
        % inter-bout interval pre and post and then get the spike trains
        % Also, get the spontaneous window spike count, the pre-bout window
        % spike count, the PST and the raster for each trial. Also, get the
        % identities of the first syllable. 
        % This can all be saved to a .mat file for faster loading and for
        % all further analysis
        disp('Calculating spike counts, rasters and psts for each neuron ...');
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.', SURecordingDetails(i).Neurontype, '.InterboutInterval', num2str(InterBoutInterval), '.SpikeChanNo', num2str(SURecordingDetails(i).SpikeChanNo), '.BoutSpikingDetails.mat']);
        if (exist(SavedFileName, 'file'))
            load(SavedFileName);
            SURecordingDetails(i).UnDirBoutSpikingDetails = Temp(1).UnDirBoutSpikingDetails;
            SURecordingDetails(i).DirBoutSpikingDetails = Temp(1).DirBoutSpikingDetails;
        else
            [SURecordingDetails(i).UnDirBoutSpikingDetails, SURecordingDetails(i).DirBoutSpikingDetails] = LSINA_CalcSpikeCountRasterPST(SURecordingDetails(i));        
            clear Temp;
            Temp(1).UnDirBoutSpikingDetails = SURecordingDetails(i).UnDirBoutSpikingDetails;
            Temp(1).DirBoutSpikingDetails = SURecordingDetails(i).DirBoutSpikingDetails;
            save(SavedFileName, 'Temp');
        end

        % Now for each neuron, find all the bouts that have the required
        % inter-bout interval pre and post and then get the spike trains
        % Also, get the spontaneous window spike count, the post-bout window
        % spike count, the PST and the raster for each trial. Also, get the
        % identities of the first syllable. 
        % This can all be saved to a .mat file for faster loading and for
        % all further analysis
        % This is all for bout offset, while teh previous set for bout
        % onset
        disp('Calculating spike counts, rasters and psts for each neuron ...');
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.', SURecordingDetails(i).Neurontype, '.InterboutInterval', num2str(InterBoutInterval), '.SpikeChanNo', num2str(SURecordingDetails(i).SpikeChanNo), '.BoutSpikingDetails_BoutOffset.mat']);
        if (exist(SavedFileName, 'file'))
            load(SavedFileName);
            SURecordingDetails(i).UnDirBoutSpikingDetails_BoutOffset = Temp(1).UnDirBoutSpikingDetails_BoutOffset;
            SURecordingDetails(i).DirBoutSpikingDetails_BoutOffset = Temp(1).DirBoutSpikingDetails_BoutOffset;
        else
            [SURecordingDetails(i).UnDirBoutSpikingDetails_BoutOffset, SURecordingDetails(i).DirBoutSpikingDetails_BoutOffset] = LSINA_CalcSpikeCountRasterPST_BoutOffset(SURecordingDetails(i));        
            clear Temp;
            Temp(1).UnDirBoutSpikingDetails_BoutOffset = SURecordingDetails(i).UnDirBoutSpikingDetails_BoutOffset;
            Temp(1).DirBoutSpikingDetails_BoutOffset = SURecordingDetails(i).DirBoutSpikingDetails_BoutOffset;
            save(SavedFileName, 'Temp');
        end
        
        % Now calculate gap spontaneous activity, basically find all
        % silent gaps that are 0.5s long and have >= 3s of silence
        % before and after. This would then be used as a measure of
        % spontaneous activity.
        disp('Calculating gap spontaneous activity for each neuron ...');
        SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.', SURecordingDetails(i).Neurontype, '.InterboutInterval', num2str(InterBoutInterval), '.SpikeChanNo', num2str(SURecordingDetails(i).SpikeChanNo), '.GapSpikingDetails.mat']);
        if (exist(SavedFileName, 'file'))
            load(SavedFileName);
            SURecordingDetails(i).UnDirGapSpontActivity = Temp(1).UnDirGapSpontActivity;
            SURecordingDetails(i).DirGapSpontActivity = Temp(1).DirGapSpontActivity;
        else
            [SURecordingDetails(i).UnDirGapSpontActivity, SURecordingDetails(i).DirGapSpontActivity] = LSINA_CalcSilentGapSpontActivity(SURecordingDetails(i));
            clear Temp;
            Temp(1).UnDirGapSpontActivity = SURecordingDetails(i).UnDirGapSpontActivity;
            Temp(1).DirGapSpontActivity = SURecordingDetails(i).DirGapSpontActivity;
            save(SavedFileName, 'Temp');
        end
        
        disp('Putting together list of syllables ...');
        Saved_FileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.', SURecordingDetails(i).Neurontype, '.InterboutInterval', num2str(InterBoutInterval), '.SpikeChanNo', num2str(SURecordingDetails(i).SpikeChanNo), '.SyllableInfo.mat']);
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
        Saved_FileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.', SURecordingDetails(i).Neurontype, '.InterboutInterval', num2str(InterBoutInterval), '.SpikeChanNo', num2str(SURecordingDetails(i).SpikeChanNo), '.SAPFeats.mat']);
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
            SAPFeatsMatrix = SURecordingDetails(i).SAPFeatsMatrix;
            SAPFeat_FieldNames = SURecordingDetails(i).SAPFeat_FieldNames;
            save(Saved_FileName, 'SAPFeatsMatrix', 'SAPFeat_FieldNames');
        end

        % Now to get an index that gives the cumulative syll # for each of
        % the bouts
        disp('Calculating cumulative syllable # for each bout ...');
        Saved_FileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.', SURecordingDetails(i).Neurontype, '.InterboutInterval', num2str(InterBoutInterval), '.SpikeChanNo', num2str(SURecordingDetails(i).SpikeChanNo), '.CumulSyllNumBoutWise.mat']);
        if (exist(Saved_FileName, 'file'))
            load(Saved_FileName);
            SURecordingDetails(i).BoutSyllNo = BoutSyllNo;
        else
            CumulSyllNo = 0;
            for j = 1:size(SURecordingDetails(i).Bouts,1),
                SURecordingDetails(i).BoutSyllNo(j,:) = CumulSyllNo + [1 (SURecordingDetails(i).Bouts(j,2) - SURecordingDetails(i).Bouts(j,1) + 1)];
                CumulSyllNo = SURecordingDetails(i).BoutSyllNo(j,end);
            end
            BoutSyllNo = SURecordingDetails(i).BoutSyllNo;
            save(Saved_FileName, 'BoutSyllNo');
        end
        
        % Now calculate IFR (instantaneous firing rates) for each valid
        % bout and save that. Do this for hte spontaneous activity gaps too
        
%         disp('Calculating IFR for each neuron ...');
%         SavedFileName = fullfile(SavedFileDir, [SURecordingDetails(i).SongFileList, '.', SURecordingDetails(i).Neurontype, '.InterboutInterval', num2str(InterBoutInterval), '.IFR.mat']);
%         if (exist(SavedFileName, 'file'))
%             load(SavedFileName);
%             SURecordingDetails(i).UnDirBoutSpikingDetailsIFR = Temp(1).UnDirBoutSpikingDetailsIFR;
%             SURecordingDetails(i).DirBoutSpikingDetailsIFR = Temp(1).DirBoutSpikingDetailsIFR;
%             SURecordingDetails(i).UnDirGapSpontActivityIFR = Temp(1).UnDirGapSpontActivityIFR;
%             SURecordingDetails(i).DirGapSpontActivityIFR = Temp(1).DirGapSpontActivityIFR;
%         else
%             [SURecordingDetails(i).UnDirBoutSpikingDetailsIFR, SURecordingDetails(i).DirBoutSpikingDetailsIFR, SURecordingDetails(i).UnDirGapSpontActivityIFR, SURecordingDetails(i).DirGapSpontActivityIFR] = LSINA_CalcIFR(SURecordingDetails(i));
%             clear Temp;
%             Temp(1).UnDirBoutSpikingDetailsIFR = SURecordingDetails(i).UnDirBoutSpikingDetailsIFR;
%             Temp(1).DirBoutSpikingDetailsIFR = SURecordingDetails(i).DirBoutSpikingDetailsIFRTemp(1).DirBoutSpikingDetailsIFR;
%             Temp(1).UnDirGapSpontActivityIFR = SURecordingDetails(i).UnDirGapSpontActivityIFR;
%             Temp(1).DirGapSpontActivityIFR = SURecordingDetails(i).DirGapSpontActivityIFR;
%             save(SavedFileName, 'Temp');
%         end
        
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
    
%     disp('Calculating IFRs ...');
%     % I am calculating the IFR for each trial and then convolving it with a
%     % gaussian of 5ms width
%     GaussianWidth = 5; % ms
%     GaussianLen = 1; % in std
%     IFR_Fs = 1000; % IFR sampling freq is 1000 Hz
%     
%     for i = 1:length(SURecordingDetails),
%         [SURecordingDetails(i).ClusterIFR, SURecordingDetails(i).OptionalClusterIFR] = LSINA_CalculateSmoothedIFRs(SURecordingDetails(i), GaussianWidth/1000, GaussianLen, IFR_Fs);
%     end
%     
%     disp('Splitting files into directed and undirected...');
%     for i = 1:length(SURecordingDetails),
%         [SURecordingDetails(i).SpikeInfo] = LSINA_LoadSpikeTimes(SURecordingDetails(i));
%     end
%     
% end
