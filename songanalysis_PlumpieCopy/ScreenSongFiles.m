function [ActualSongFiles, NoSongFiles, ProblemSongFiles] = ScreenSongFiles(DirectoryName,BirdName,FileType)

if (DirectoryName(end) ~= '/')
    DirectoryName(end + 1) = '/';
end

cd(DirectoryName);

AllSongFiles = dir(['*',BirdName,'*']);

ActualSongFiles = [];
AIndex = 0;
ProblemSongFiles = [];
PIndex = 0;
NoSongFiles = [];
NIndex = 0;

for i = 1:length(AllSongFiles),
    %figure;
    try
        PlotSpectrogram(DirectoryName,AllSongFiles(i).name,FileType, 'hot');
        zoom xon;
        SongOrNoSong = menu(['File #',num2str(i),' of ',num2str(length(AllSongFiles)),': Has Song or Does Not Have Song'],'File has song','File does not have song', 'Play File', 'Quit');
        if (SongOrNoSong == 1)
            AIndex = AIndex + 1;
            ActualSongFiles{AIndex} = [DirectoryName, AllSongFiles(i).name];
            %fid = fopen('ActualSongFileNames.txt','a');
            for j = 1:length(ActualSongFiles{AIndex});
                %fprintf(fid,'%c',ActualSongFiles{AIndex}(j));
            end
            %fprintf(fid,'\n');
            %fclose(fid);
        else
            if (SongOrNoSong == 4)
                 c;
                 break;
            else
                if (SongOrNoSong == 3)
                    [rawsong,Fs] = soundin_copy(DirectoryName, AllSongFiles(i).name,'obs0r');
                    %rawsong = rawsong/32768;
                    soundsc(rawsong, Fs);
                    i = i - 1;
                else
                    NIndex = NIndex + 1;
                    NoSongFiles{NIndex} = AllSongFiles(i).name;
                end
            end
        end
    catch
        PIndex = PIndex + 1;
        ProblemSongFiles{PIndex} = AllSongFiles(i).name;
        %fid = fopen('ProblemSongFileNames.txt','a');
        for j = 1:length(ProblemSongFiles{PIndex});
            %fprintf(fid,'%c',ProblemSongFiles{PIndex}(j));
        end
        %fprintf(fid,'\n');
        %fclose(fid);
    end
    close(gcf);
end

