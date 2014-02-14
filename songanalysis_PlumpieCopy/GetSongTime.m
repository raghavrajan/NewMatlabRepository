function [FileTime] = GetSongTime(FileName)


Cbin_index = strfind(FileName,'.cbin');
if (length(Cbin_index) > 0)
    RecFilename = [FileName(1:Cbin_index),'rec'];
    fid2 = fopen(RecFilename);
    while (~(feof(fid2)))
        tline2 = fgetl(fid2);
        if (strfind(tline2,'File created'))
            Commaindex = strfind(tline2,',');
            Date = tline2((Commaindex(1) - 3):(Commaindex(end) - 1));
            Time = tline2((Commaindex(end) + 1):end);
            Colonindex = strfind(Time,':');
            FileTime = str2double(Time(1:(Colonindex(1) - 1))) + (str2double(Time((Colonindex(1)+1):(Colonindex(2) - 1))))/60 + (str2double(Time((Colonindex(2) + 1):end))/3600);
        end
        if (strfind(tline2,'rec end'))
            EqualIndex = strfind(tline2,'=');
            RecordLength = str2double(tline2((EqualIndex + 1):(end - 2)))/1000;
        end
    end
    fclose(fid2);
else
    FileTime = 0;
%     DotIndex = strfind(FileName,'.');
%     fid2 = fopen([FileName(1:(DotIndex(2) - 1)),'.log']);
%     while (~feof(fid2))
%        tline = fgetl(fid2);
%        strfind(tline,FileName(1:(DotIndex(3) - 1)));
%        
%     end
%     fclose(fid2);
end
