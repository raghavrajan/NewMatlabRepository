function [] = InteractiveDoAnalysis()

RawDataDirectory = uigetdir(pwd, 'Choose the data directory');
cd(RawDataDirectory);
GenParam.RawDataDirectory = RawDataDirectory;

RootFileName = inputdlg('Enter the root file name for all the data files');
RootFileName = RootFileName{1};

[DirFileNames] = uigetfile({[RootFileName,'*'], 'All Files (*.*)'}, 'Choose the directed song files', 'MultiSelect', 'on');
[UnDirFileNames] = uigetfile({[RootFileName,'*'], 'All Files (*.*)'}, 'Choose the undirected song files', 'MultiSelect', 'on');

if (length(DirFileNames) == 1)
    DirFileInfo.FileNames{1} = DirFileNames;
else
    DirFileInfo.FileNames = DirFileNames;
end
clear DirFileNames;

if (length(UnDirFileNames) == 1)
    UnDirFileInfo.FileNames{1} = UnDirFileNames;
else
    UnDirFileInfo.FileNames = UnDirFileNames;
end
clear UnDirFileNames;

FileType = inputdlg('Enter the file type - obs for Observer files, okrank for OKrank files and wav for Wav files');
FileType = FileType{1};
GenParam.FileType = FileType;

[DirFileInfo.RecordLengths, DirFileInfo.NChans] = GetRecordLengths(RawDataDirectory, DirFileInfo.FileNames, FileType);
[UnDirFileInfo.RecordLengths, UnDirFileInfo.NChans] = GetRecordLengths(RawDataDirectory, UnDirFileInfo.FileNames, FileType);

NotesFilesDirectory = uigetdir(pwd, 'Choose the notes files directory');
GenParam.NotesFilesDirectory = NotesFilesDirectory;

[DirFileInfo.NoteOnsets, DirFileInfo.NoteOffsets, DirFileInfo.NoteLabels] = GetNoteInformation(GenParam.NotesFilesDirectory, DirFileInfo.FileNames);
[UnDirFileInfo.NoteOnsets, UnDirFileInfo.NoteOffsets, UnDirFileInfo.NoteLabels] = GetNoteInformation(GenParam.NotesFilesDirectory, UnDirFileInfo.FileNames);

Motif = inputdlg('Enter the motif that has to be analysed');
Motif = Motif{1};
GenParam.Motif = Motif;

[DirFileInfo.Syllables, DirFileInfo.Gaps] = CalculateSyllableStatistics(DirFileInfo, GenParam.Motif);
[UnDirFileInfo.Syllables, UnDirFileInfo.Gaps] = CalculateSyllableStatistics(UnDirFileInfo, GenParam.Motif);

SpikeSortMethod = inputdlg('Enter the method to be used for spike sorting - 0 for thresholding with positive first and then negative for spikes, 1 for thresholding with negative first and then positive and 2 for already sorted files');
SpikeSortMethod = SpikeSortMethod{1};

GenParam.SpikeSortMethod = SpikeSortMethod;

SpikeChanNo = inputdlg('Enter the channel no for the spike data - 0 being the first channel');
SpikeChanNo = SpikeChanNo{1};
GenParam.SpikeChanNo = SpikeChanNo;

if (SpikeSortMethod < 2)
    [DirFileInfo.SpikeData, GenParam.ThresholdingParameters] = GetSpikeTimes(GenParam, DirFileInfo);    
    [UnDirFileInfo.SpikeData, GenParam.ThresholdingParameters] = GetSpikeTimes(GenParam, UnDirFileInfo);    
end
    
disp('Finished');
