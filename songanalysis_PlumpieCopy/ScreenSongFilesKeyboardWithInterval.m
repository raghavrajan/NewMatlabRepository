function [ActualSongFiles, ActualCallFiles, BaselineFiles, ProblemSongFiles] = ScreenSongFilesKeyboardWithInterval(DirectoryName, BirdName, FileType, OutputFileName, InterBoutInterval)

if (DirectoryName(end) ~= '/')
    DirectoryName(end + 1) = '/';
end

cd(DirectoryName);

AllSongFiles = dir([BirdName,'*']);

ActualSongFiles = [];
AIndex = 0;
ActualCallFiles = [];
CIndex = 0;
ProblemSongFiles = [];
PIndex = 0;
BaselineFiles = [];
BIndex = 0;

for i = 1:length(AllSongFiles),
    %figure('Position',[6 451 1270 317]);
    try
        PlotLowResSpectrogram(DirectoryName,AllSongFiles(i).name,FileType);
        hold on;
        plot([0 InterBoutInterval], [7000 7000], 'b', 'LineWidth', 2);
        Temp = axis;
        plot([(Temp(2) - InterBoutInterval) Temp(2)], [7000 7000], 'b', 'LineWidth', 2);
        
        title([AllSongFiles(i).name, ' #', num2str(i), ' of ', num2str(length(AllSongFiles)), ' files: SongFile # ', num2str(AIndex)]);
        zoom xon;
        [x, y, SongOrNoSong] = ginput(1);
        if (SongOrNoSong == 115)
            AIndex = AIndex + 1;
            ActualSongFiles{AIndex} = AllSongFiles(i).name;
        else
            if (SongOrNoSong == 113)
                 break;
            else
                if (SongOrNoSong == 99)
                    CIndex = CIndex + 1;
                    ActualCallFiles{CIndex} = AllSongFiles(i).name;
                else
                    if (SongOrNoSong == 98)
                        BIndex = BIndex + 1;
                        BaselineFiles{BIndex} = AllSongFiles(i).name;
                    end
                end
            end
        end
    catch
        PIndex = PIndex + 1;
        ProblemSongFiles{PIndex} = AllSongFiles(i).name;
        fid = fopen('ProblemSongFileNames.txt','a');
        for j = 1:length(ProblemSongFiles{PIndex});
            fprintf(fid,'%c',ProblemSongFiles{PIndex}(j));
        end
        fprintf(fid,'\n');
        fclose(fid);
    end
    close(gcf);
end

Fid = fopen(OutputFileName, 'w');
for i = 1:length(ActualSongFiles),
    fprintf(Fid, '%s\n', ActualSongFiles{i});
end
fclose(Fid);