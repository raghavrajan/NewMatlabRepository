function [] = FindAndWriteSongFilesToFileList(BirdDetailsTextFile, InterboutInterval)

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

% Now based on the bouts identified, write out the list of files with song
% bouts in them with enough time - if continuous data, then just write
% everything out
for i = 1:length(BirdParameters),
    OutputFileName = fullfile(BirdParameters(i).DataDirectory, BirdParameters(i).SongFileList, '.SongFiles.txt');
    if (exist(OutputFileName, 'file'))
        disp('Output song file list already exists');
        continue;
    else
        Fid = fopen(fullfile(BirdParameters(i).DataDirectory, [BirdParameters(i).SongFileList, '.SongFiles.txt']), 'w');
        ValidBouts = find(BirdParameters(i).Bouts(:,7) == 1);
        if (BirdParameters(i).Continuousdata == 1)
            ValidFileNumbers = [BirdParameters(i).Bouts(ValidBouts,3); BirdParameters(i).Bouts(ValidBouts,4)];
        else
            ValidFileNumbers = BirdParameters(i).Bouts(ValidBouts,3);
        end
        ValidFileNumbers = unique(ValidFileNumbers);
        for j = 1:length(ValidFileNumbers),
            fprintf(Fid, '%s\n', BirdParameters(i).SongFileNames{ValidFileNumbers(j)});
        end
        fclose(Fid);
    end
end
disp('Finished writing song file lists');



