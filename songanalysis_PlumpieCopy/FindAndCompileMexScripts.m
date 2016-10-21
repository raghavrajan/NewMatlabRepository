function [] = FindAndCompileCFiles

cd ~/repositories/Rajan-Lab-Song-Analysis-Scripts/songanalysis_PlumpieCopy

Files = dir('*.c');
Files = dir('*');
for i = 1:length(Files),
    
    if isdir(Files(i).name);
    try
        eval(['mex ', Files(i).name]);
    catch
        disp(['Could not compile ', Files(i).name]);
    end
end
