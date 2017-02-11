function [ActualDataFiles, ProblemDataFiles] = ScreenDataFiles(DirectoryName,BirdName,FileType)

if (DirectoryName(end) ~= '/')
    DirectoryName(end + 1) = '/';
end

cd(DirectoryName);

AllDataFiles = dir([BirdName,'*']);

ActualDataFiles = [];
AIndex = 0;
ProblemDataFiles = [];
PIndex = 0;
for i = 1:length(AllDataFiles),
    %figure('Position',[6 451 1270 317]);
    try
        if (strfind(AllDataFiles(i).name, '.rec'))
            continue;
        end
        PlotOKrankData(DirectoryName,AllDataFiles(i).name);
        title(AllDataFiles(i).name);
        zoom xon;
        DataOrNoData = menu(['File #',num2str(i),' of ',num2str(length(AllDataFiles)),': Has Data or Does Not Have Data'],'File has Data','File does not have Data', 'Quit');
        if (DataOrNoData == 1)
            AIndex = AIndex + 1;
            ActualDataFiles{AIndex} = AllDataFiles(i).name;
            fid = fopen('ActualDataFileNames.txt','a');
            for j = 1:length(ActualDataFiles{AIndex});
                fprintf(fid,'%c',ActualDataFiles{AIndex}(j));
            end
            fprintf(fid,'\n');
            fclose(fid);
        else
             if (DataOrNoData == 3)
                 c;
                 break;
             end
        end
    catch
        PIndex = PIndex + 1;
        ProblemDataFiles{PIndex} = AllDataFiles(i).name;
        fid = fopen('ProblemDataFileNames.txt','a');
        for j = 1:length(ProblemDataFiles{PIndex});
            fprintf(fid,'%c',ProblemDataFiles{PIndex}(j));
        end
        fprintf(fid,'\n');
        fclose(fid);
    end
    close(gcf);
end

