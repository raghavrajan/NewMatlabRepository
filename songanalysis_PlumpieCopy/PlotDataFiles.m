function [] = PlotDataFiles(DirectoryName, FileList, ChanNo)

if ~(exist('FileList1','var'))
    fid = fopen(FileList, 'r');
    FileName = fgetl(fid);

    while (ischar(FileName(1)))
        PlotOKrankData(DirectoryName, FileName, ChanNo);
        FileName = fgetl(fid);
    end
    fclose(fid);
else
    for i = 1:length(FileList1),
        PlotOKrankData(DirectoryName, FileList1{i}, ChanNo);
    end
end
