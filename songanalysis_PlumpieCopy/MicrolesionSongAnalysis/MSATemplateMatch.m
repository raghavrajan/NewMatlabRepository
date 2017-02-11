function [] = MSATemplateMatch(RawDataDir, SongFileList, FileType, NotesFile)

PresentDir = pwd;
Notes = load(NotesFile);
SongFile = NotesFile(1:end-8);

[Template, TemplateLabels] = MSAMakeAllSyllableTemplatesFromFile(RawDataDir, PresentDir, SongFile, FileType, PresentDir, 0);

disp('Finished template matching');